from pathlib import Path
import os
import logging
import sys

logger = logging.getLogger(__name__)

from gdex_api_client import gdex_client as gac

ENVAUTHFILE = 'GLOBSIM_RDA_AUTH_FILE'

def lookup_auth_file(explicit_path=None):
    """Look for an authentication file in order of precedence:
    1. Explicit path passed to method/class
    2. Environment variable GLOBSIM_RDA_AUTH_FILE
    3. Default location: ~/rdams_token.txt
    """
    if explicit_path:
        p = Path(explicit_path).expanduser()
        if p.exists():
            return str(p)

    env_auth = os.environ.get(ENVAUTHFILE, None)
    if env_auth and Path(env_auth).expanduser().exists():
        return str(Path(env_auth).expanduser())

    default_auth = Path("~", 'rdams_token.txt').expanduser()
    if default_auth.exists():
        return str(default_auth)

    return None


def globsim_get_authentication(auth_file=None):
    """Get the authentication token using the resolved auth file."""
    resolved_file = lookup_auth_file(auth_file)
    if not resolved_file:
        logger.error(
            f"No authentication file found. Please set {ENVAUTHFILE} "
            "or place a token file at ~/rdams_token.txt."
        )
        sys.exit(1)

    return gac.read_token_file(resolved_file)


# Patch gdex_client directly
gac.get_authentication = globsim_get_authentication


import os
import sys
import time
import logging
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

logger = logging.getLogger(__name__)

def fast_download_files(filelist, out_dir='./', cookie_file=None, max_workers=2, max_retries=3):
    """Robust replacement for gac.download_files with connection pooling,
    browser header spoofing, and automatic socket-reset on throttle.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:120.0) Gecko/20100101 Firefox/120.0',
        'Accept': '*/*',
        'Connection': 'keep-alive'
    })

    def download_single(url):
        filename = url.split('/')[-1]
        dest = out_path / filename
        temp_dest = out_path / f"{filename}.part"

        # 1. Skip if the final completed file already exists
        if dest.exists() and dest.stat().st_size > 0:
            logger.info(f"Skipping existing file: {filename}")
            return

        for attempt in range(1, max_retries + 1):
            try:
                print(f"Downloading {filename} (Attempt {attempt})...")
                
                # Download to temporary file .part
                with session.get(url, stream=True, timeout=(10, 15)) as r:
                    r.raise_for_status()
                    with open(temp_dest, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                f.write(chunk)

                # 2. Rename .part to final destination ONLY on 100% completion
                temp_dest.replace(dest)
                print(f"100% Completed: {filename}")
                return

            except (requests.exceptions.Timeout, requests.exceptions.RequestException) as e:
                logger.warning(f"Connection stalled on {filename} ({e}). Cleaning up temp file...")
                
                # Clean up the partial .part file if download failed mid-way
                if temp_dest.exists():
                    temp_dest.unlink()
                    
                time.sleep(2)

        print(f"Error: Failed to download {filename} after {max_retries} attempts.")

    # Execute downloads across thread pool
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(download_single, set(filelist)))


# Patch the API package function directly
gac.download_files = fast_download_files


class MetaRdams(type):
    """Metaclass allowing class-level attribute lookup delegation to gac."""
    def __getattr__(cls, name):
        try:
            return getattr(gac, name)
        except AttributeError:
            raise AttributeError(f"Neither 'Rdams' nor 'gdex_client' has attribute '{name}'")


class Rdams(metaclass=MetaRdams):
    """Thin wrapper around gdex_client (gac). 
    Delegates all undefined class or instance attributes/methods directly to gac.
    """
    ENVAUTHFILE = ENVAUTHFILE

    def __init__(self, auth_file=None):
        self.auth_file = auth_file
        if auth_file:
            path = Path(auth_file).expanduser()
            if path.exists():
                os.environ[self.ENVAUTHFILE] = str(path)
            else:
                logger.warning(f"Provided auth file does not exist: {auth_file}")

    def __getattr__(self, name):
        """Instance-level delegation to gac."""
        try:
            return getattr(gac, name)
        except AttributeError:
            raise AttributeError(f"Neither 'Rdams' nor 'gdex_client' has attribute '{name}'")

def parse_rinfo(rinfo):
    """Parse the rinfo string into a dictionary of parameters."""
    param_dict = {}
    for item in rinfo.split(';'):
        key, value = item.split('=')
        param_dict[key.strip()] = value.strip()
    return param_dict

def parse_rinfo_parameters(param_string):
    """ '8!d640000:tprate1have-sfc-fc-gauss,8!d640000:dswrf1have-sfc-fc-gauss,8!d640000:dlwrf1have-sfc-fc-gauss,8!d640000:dswrfcs1have-sfc-fc-gauss,8!d640000:dlwrfcs1have-sfc-fc-gauss'
    Parse the parameter string into a list.
    """
    param_list = param_string.split(',')
    parsed_params = []
    for param in param_list:
        if ':' in param:
            dataset, variable = param.split(':', 1)
            parsed_params.append({'dataset': dataset, 'variable': variable})
        else:
            logger.warning(f"Unexpected parameter format: {param}")
    return parsed_params

def get_parsed_status():
    status = gac.get_status()

    if status['http_response'] != 200:
        logger.error(f"Failed to get status: {status}")
        return None
    if status.get('data') is None:
        logger.warning("No data in status response.")
        return None

    if len(data := status['data']) > 0:
        for request in data:
            rinfo = request.get('rinfo')
            if rinfo:
                parsed_info = parse_rinfo(rinfo)
                parsed_info['parsed_parameters'] = parse_rinfo_parameters(parsed_info.get('parameters', ''))
                request['parsed_rinfo'] = parsed_info
            else:
                logger.warning(f"No rinfo found for request: {request}")
    
    return status
            
    

if __name__ == "__main__":
    # Example usage:
    rdams = Rdams(auth_file=None)
    # Delegates cleanly to gac.get_summary using globsim_get_authentication
    summary = rdams.get_summary('ds640.0')
      
    status=rdams.get_status()
    rinfo = status['data'][0]['rinfo']
    import pdb; pdb.set_trace() 
    d = get_request_dicts()
    