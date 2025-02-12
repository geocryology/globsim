#!/bin/bash

#!/bin/bash

# Default values
YEAR=""
TIME="valid_time"
LEV="pressure_level"
PRODUCT="era5"
OVERWRITE=false

# 'Hard-coded' values
LAT='latitude'
LON='longitude' 
LAT_CHK=6
LON_CHK=6

# Help message
usage() {
    echo "Usage: $0 input_dir output_dir [--year YEAR] [--time TIME] [--lev LEV] [--overwrite]"
    echo
    echo "Positional arguments:"
    echo "  input_dir   Path to input directory"
    echo "  output_dir  Path to output directory"
    echo
    echo "Keyword arguments:"
    echo "  --product PRODUCT Specify the reanalysis (era5, jra55, etc.)"
    echo "  --year YEAR       Specify the year to process (default: empty)"
    echo "  --time TIME       Variable name for time dimension (default: valid_time)"
    echo "  --lev LEV         Variable name for level dimension (default: pressure_level)"
    echo "  --overwrite       Overwrite existing files (default: false)"
    echo
    exit 1
}

# Check if at least two arguments are provided
if [ "$#" -lt 2 ]; then
    usage
fi

# Positional arguments
input_dir="$1"
output_dir="$2"
shift 2

# Parse keyword arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --year)
            YEAR="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --product)
            PRODUCT="$2"
            shift 2
            ;;
        --lev)
            LEV="$2"
            shift 2
            ;;
        --overwrite)
            OVERWRITE=true
            shift
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Print parsed values
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"
echo "Year: ${YEAR:-Not specified}"
echo "Time variable: $TIME"
echo "Level variable: $LEV"
echo "Product variable: $PRODUCT"
echo "Overwrite: $OVERWRITE"

get_dim_size() {
    ncdump -h "$2" | awk -v dim="$1" '$1 == dim && $2 == "=" {print $3}' | tr -d ';'
}

# Define dimension orders for each type
declare -A dim_orders

dim_orders[pl]="latitude,longitude,pressure_level,valid_time"
dim_orders[sa]="latitude,longitude,valid_time"
dim_orders[sf]="latitude,longitude,valid_time"

# ensure paths are different
if [ "$(realpath "$input_dir")" == "$(realpath "$output_dir")" ]; then
  echo "Input and output paths are identical"
  exit 1
fi

# Make sure output directory exists
mkdir -p "$output_dir"

# Loop over all NetCDF files in the input directory
for file in "$input_dir"/${PRODUCT}_*_${YEAR}*_to_*.nc; do
    
    echo $file
    # Extract the filename from the full path
    filename=$(basename "$file")
    
    TIME_CHK=$(get_dim_size $TIME $file)
    LEV_CHK=$(get_dim_size $LEV $file)
    # Extract the type of file (pl, sa, sf) from the filename
    if [[ $filename =~ ${PRODUCT}_([a-z]+)_ ]]; then
        type=${BASH_REMATCH[1]}
        if [[ -n ${dim_orders[$type]} ]]; then
            # Define the output file path
            output_file="$output_dir/$filename"
            orig_file="$input_dir/orig_$filename"
            if [[ -f "$output_file" && "$OVERWRITE" -eq 0 ]]; then  
                echo "Skipping file: $file (OVERWRITE is 0)"
                continue  # Skip to the next iteration
            fi

            declare -A cnk_dim
              cnk_dim[pl]="--cnk_dmn=$LON,$LON_CHK --cnk_dmn=$LAT,$LAT_CHK --cnk_dmn=$TIME,$TIME_CHK --cnk_dmn=$LEV,$LEV_CHK"
              cnk_dim[sa]="--cnk_dmn=$LON,$LON_CHK --cnk_dmn=$LAT,$LAT_CHK --cnk_dmn=$TIME,$TIME_CHK"
              cnk_dim[sf]="--cnk_dmn=$LON,$LON_CHK --cnk_dmn=$LAT,$LAT_CHK --cnk_dmn=$TIME,$TIME_CHK"
            
            #echo "Reordering $file -> $output_file with order: ${dim_orders[$type]}"
            #ncpdq --rdr=${dim_orders[$type]} ${cnk_dim[$type]} "$file" "$output_file"
            time ncpdq --rdr=${dim_orders[$type]} ${cnk_dim[$type]} "$file" "$output_file"  # provide timing information
            echo "ncpdq --rdr=${dim_orders[$type]} ${cnk_dim[$type]} \"$file\" \"$output_file\""  # just print the command
        else
            echo "Skipping $filename: Unknown type"
        fi
    fi
done
