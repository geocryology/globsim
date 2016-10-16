#
# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.
#

import os
import sys
import time
import urllib
import urllib2
import httplib

try:
    import json
except:
    import simplejson as json

try:
    import socket
    socket.ssl
except:
    print "Python socket module was not compiled with SSL support. Aborting..."
    sys.exit(1)

###############################################################################
VERSION = '1.2'

###############################################################################

KEY = None
URL = None
EMAIL = None

rc = os.path.normpath(os.path.join(os.environ["HOME"],".ecmwfapirc"))
if os.path.exists(rc):
    try:
        config = json.loads(file(rc).read())
        URL   = config.get("url", "https://api.ecmwf.int/v1")
        KEY   = config.get("key")
        EMAIL = config.get("email")
    except Exception, e:
        print e
        sys.exit(1)

KEY=os.environ.get("ECMWF_API_KEY", KEY)
URL=os.environ.get("ECMWF_API_URL", URL)
EMAIL=os.environ.get("ECMWF_API_EMAIL", EMAIL)

###############################################################################

class RetryError(Exception):
    def __init__(self, code, text):
       self.code = code
       self.text = text
    def __str__(self):
        return "%d %s" % (self.code, self.text)

class APIException(Exception):
    def __init__(self, value):
       self.value = value
    def __str__(self):
        return repr(self.value)

def robust(func):

    def wrapped(*args,**kwargs):
        tries = 0
        while True:
            try:
                return func(*args,**kwargs)
            except urllib2.HTTPError, e:
                print "WARNING: httplib2.HTTPError received %s" % (e)
                if e.code < 500: raise
                tries += 1
                if tries > 10: raise
                time.sleep(60)
            except httplib.BadStatusLine, e:
                print "WARNING: httplib.BadStatusLine received %s" % (e)
                tries += 1
                if tries > 10: raise
                time.sleep(60)
            except urllib2.URLError, e:
                print "WARNING: httplib2.URLError received %s %s" % (e.errno, e)
                tries += 1
                if tries > 10: raise
                time.sleep(60)
            except APIException:
                raise
            except RetryError, e:
                print "WARNING: HTTP received %s" % (e.code)
                print e.text
                tries += 1
                if tries > 10: raise
                time.sleep(60)
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise

    return wrapped

SAY = True
class Ignore303(urllib2.HTTPRedirectHandler):


    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code in [301, 302]:
            # We want the posts to work even if we are redirected
            if code == 301:
                global SAY, URL
                if SAY:
                    o = req.get_full_url()
                    n = newurl
                    while o != URL and len(o) and len(n) and o[-1] == n[-1]:
                        o = o[0:-1]
                        n = n[0:-1]
                    print
                    print "*** ECMWF API has moved"
                    print "***   OLD: %s" % o
                    print "***   NEW: %s" % n
                    print "*** Please update your ~/.ecmwfapirc file"
                    print
                    SAY = False
            data = None
            if req.has_data():
                data = req.get_data()
            return urllib2.Request(newurl, data=data, headers = req.headers, 
                                    origin_req_host=req.get_origin_req_host(), unverifiable=True)
        return None

    def http_error_303(self, req, fp, code, msg, headers):
        infourl = urllib.addinfourl(fp, headers, req.get_full_url())
        infourl.status = code
        infourl.code = code
        return infourl

class Connection(object):

    def __init__(self, email = None, key = None, verbose = False, quiet = False):
        self.email    = email
        self.key      = key
        self.retry    = 5
        self.location = None
        self.done     = False
        self.value    = True
        self.offset   = 0
        self.verbose  = verbose
        self.quiet    = quiet
        self.status   = None

    @robust
    def _call(self, url, payload = None, method = "GET"):

        if self.verbose:
            print method, url

        headers = { "Accept" : "application/json", "From" : self.email, "X-ECMWF-KEY" : self.key }

        opener = urllib2.build_opener(Ignore303)

        data = None
        if payload is not None:
            data = json.dumps(payload)
            data.encode('utf-8')
            headers["Content-Type"] = "application/json";

        url = "%s?offset=%d&limit=500" % (url, self.offset)
        req = urllib2.Request(url=url, data=data, headers=headers)
        if method:
            req.get_method = lambda: method

        error = False
        try:
            try:
                res  = opener.open(req)
            except urllib2.HTTPError,e:
                # It seems that some version of urllib2 are buggy
                if e.code <= 299:
                    res = e
                else:
                    raise
        except urllib2.HTTPError,e:
            print e
            error = True
            res   = e
            # 502: Proxy Error
            # 503: Service Temporarily Unavailable
            if e.code >= 500:
                raise RetryError(e.code, e.read())

        self.retry    = int(res.headers.get("Retry-After", self.retry))
        code          = res.code
        if code in [201, 202]:
            self.location = res.headers.get("Location",    self.location)

        if self.verbose:
            print res.headers.get("Content-Type")
            print res.headers.get("Content-Length")

        body = res.read()
        res.close()


        if code in [204]:
            self.last = None
            return None
        else:
            try:
                self.last  =  json.loads(body)
            except Exception, e:
                self.last = { "error" : "%s: %s" % (e, body) }
                error = True

        if self.verbose:
            print json.dumps(self.last,indent=4)

        self.status = self.last.get("status", self.status)

        if "messages" in self.last:
            for n in self.last["messages"]:
                if not self.quiet:
                    print n
                self.offset += 1

        if code in [303]:
            self.value = self.last
            self.done  = True

        if "error" in self.last:
            #self.done   = True
            raise APIException("ecmwf.API error 1: %s" % (self.last["error"],) )

        if error:
            #self.done   = True
            raise APIException("ecmwf.API error 2: %s" % (res, ) )

        return self.last

    def submit(self, url, payload):
        self._call(url, payload, "POST")

    def wait(self):
        if self.verbose:
            print "Sleeping %s second(s)" % (self.retry)
        time.sleep(self.retry)
        self._call(self.location, None, "GET")

    def ready(self):
        return self.done

    def result(self):
        return self.value


    def cleanup(self):
        try:
            self._call(self.location, None, "DELETE")
        except:
            pass

def no_log(msg):
    pass

class Request(object):

    def __init__(self, url, service, email  = None, key  = None, log = no_log, quiet = False, verbose = False):
        self.url        = url
        self.service    = service
        self.connection = Connection(email, key, quiet = quiet, verbose = verbose)
        self.log        = log
        self.quiet      = quiet
        self.log("ECMWF API python library %s" % (VERSION,))
        self.log("ECMWF API at %s" % (self.url,))
        user = self.connection._call("%s/%s" % (self.url, "who-am-i"))
        self.log("Welcome %s" % (user["full_name"] or "user '%s'" % user["uid"],))
        try:
            news = self.connection._call("%s/%s/%s" % (self.url, self.service, "news"))
            for n in news["news"].split("\n"):
                self.log(n)
                
        except:
            pass

    def _bytename(self,size):   
        prefix = {'':'K','K':'M','M':'G','G':'T','T':'P','P':'E'}
        l    = ''
        size = size*1.0
        while 1024 < size:
            l = prefix[l]
            size = size / 1024
        s = ""
        if size > 1:
            s = "s"
        return "%g %sbyte%s" % (size,l,s)

    @robust
    def _transfer(self, url, path, size):
        self.log("Transfering %s into %s" % (self._bytename(size), path))
        self.log("From %s" % (url, ))
        start = time.time()
        http = urllib2.urlopen(url)
        f = open(path,"wb")
        total = 0
        block = 1024*1024
        while True:
            chunk = http.read(block)
            if not chunk: break
            f.write(chunk)
            total += len(chunk)
        f.flush()
        f.close()
        end = time.time()

        
        header = http.info()
        assert total == size
        length = header.get("content-length")
        if length is None:
            self.log("Warning: Content-Length missing from HTTP header")
        else:
            assert total == long(length)

        if end > start:
           self.log("Transfer rate %s/s" % self._bytename(total / ( end - start)), )

        return total


    def execute(self, request, target = None):

        status = None

        self.connection.submit("%s/%s/requests" % (self.url, self.service), request)
        if self.connection.status != status:
            status = self.connection.status
            self.log("Request is %s" % (status, ))

        while not self.connection.ready():
            if self.connection.status != status:
                status = self.connection.status
                self.log("Request is %s" % (status, ))
            self.connection.wait()

        if self.connection.status != status:
            status = self.connection.status
            self.log("Request is %s" % (status, ))

        result = self.connection.result()
        if target:
            size = -1
            tries = 0
            while size != result["size"] and tries < 10:
                size = self._transfer(result["href"], target, result["size"])
                if size != result["size"] and tries < 10:
                    tries += 1
                    self.log("Transfer interrupted, retrying...")
                    time.sleep(60)
                else:
                    break

            assert size == result["size"]

        self.connection.cleanup()

        return result
        

###############################################################################

class ECMWFDataServer(object):

    def __init__(self, url = URL, key = KEY, email = EMAIL, verbose = False, log = None):
        self.url     = url
        self.key     = key
        self.email   = email
        self.verbose = verbose
        self.log     = log

    def trace(self, m):
        if self.log:
            self.log(m)
        else:
            t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            print "%s %s" % (t,m,)

    def retrieve(self, req):
        target  = req.get("target")
        dataset = req.get("dataset")
        c = Request(URL, "datasets/%s" % (dataset,), self.email, self.key, self.trace, verbose = self.verbose)
        c.execute(req, target)

###############################################################################

class ECMWFService(object):

    def __init__(self, service, url = URL, key = KEY, email = EMAIL, verbose = False, log = None, quiet = False):
        self.service = service
        self.url     = url
        self.key     = key
        self.email   = email
        self.verbose = verbose
        self.quiet   = quiet 
        self.log     = log

    def trace(self, m):
        if self.log:
            self.log(m)
        else:
            t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            print "%s %s" % (t,m,)

    def execute(self, req, target):
        c = Request(URL, "services/%s" % (self.service,), self.email, self.key, self.trace, verbose = self.verbose, quiet = self.quiet)
        c.execute(req, target)
        self.trace("Done.")

###############################################################################
