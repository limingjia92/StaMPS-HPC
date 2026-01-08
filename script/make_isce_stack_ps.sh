#!/bin/bash

# ==============================================================================
# Script Name: make_isce_stack_ps.sh
# Description: Prepares Single Master Stack (PS) from ISCE SLCs for StaMPS.
#              Replaces 'make_single_reference_stack_isce'.
#              Process after 'stackSentinel.py'
#
# Original Logic: David Bekaert (2017)
# Refactoring & Optimization: Mingjia Li (2025)
# Language: Bash
# ==============================================================================

# --- Function: Print Usage and Exit ---
print_usage() {
    echo "Usage: $(basename "$0") reference_date [slc_stack_path] [stack_geom_path] [slc_stack_baseline_path] [rg_looks] [az_looks]"
    echo ""
    echo "Required Arguments:"
    echo "  reference_date          : Date of the single master image (YYYYMMDD)"
    echo ""
    echo "Optional Path Arguments (Defaults assume 'merged' structure):"
    echo "  slc_stack_path          : Default: 'merged/SLC'"
    echo "  stack_geom_path         : Default: 'merged/geom_reference'"
    echo "  slc_stack_baseline_path : Default: 'merged/baselines'"
    echo ""
    echo "Optional Visualization Arguments:"
    echo "  rg_looks                : Range looks for visualization (Default: 20)"
    echo "  az_looks                : Azimuth looks for visualization (Default: 5)"
    echo ""
    exit 1
}

# --- Check Required Arguments ---
if [ $# -lt 1 ]; then
    print_usage
fi

# --- Initialize Variables with Defaults ---
REF_DATE="$1"

# Handle Optional Paths with Default Values
# Usage: ${PARAMETER:-DEFAULT}
ARG_SLC="${2:-merged/SLC}"
ARG_GEOM="${3:-merged/geom_reference}"
ARG_BASE="${4:-merged/baselines}"

# Convert to Absolute Paths (readlink -f handles relative paths correctly)
SLC_STACK_PATH=$(readlink -f "$ARG_SLC")
GEOM_PATH=$(readlink -f "$ARG_GEOM")
BASELINE_PATH=$(readlink -f "$ARG_BASE")

# Visualization Parameters (Defaults for Sentinel-1 usually 20:5)
RG_LOOKS=${5:-20}
AZ_LOOKS=${6:-5}

# --- Hardcoded Parameters ---
SLC_SUFFIX=".full"
GEOM_SUFFIX=".full"
# Sentinel-1 C-band Wavelength (m)
# Note: TerraSAR-X is ~0.031, ALOS-2 is ~0.236. Change this if using other sensors.
LAMBDA=0.0554658 

# Constants
WORK_DIR=$(pwd)
INSAR_DIR="INSAR_${REF_DATE}"
GENERATE_IFGS="y"

echo "================================================================="
echo "   ISCE to StaMPS PS Stack Preparator "
echo "================================================================="
echo "Ref Date    : $REF_DATE"
echo "SLC Path    : $SLC_STACK_PATH"
echo "Geom Path   : $GEOM_PATH"
echo "Baselines   : $BASELINE_PATH"
echo "Wavelength  : $LAMBDA"
echo "Quicklooks  : Rg:$RG_LOOKS / Az:$AZ_LOOKS"
echo "================================================================="

# --- Check Reference Date Format ---
if [[ ! "$REF_DATE" =~ ^[0-9]{8}$ ]]; then
    echo "Error: reference_date must be in YYYYMMDD format."
    exit 1
fi

# --- Check if Input Directories Exist ---
if [ ! -d "$SLC_STACK_PATH" ]; then
    echo "Error: SLC path '$SLC_STACK_PATH' does not exist."
    exit 1
fi

# --- Step 1: Create Structure & Reference Setup ---
echo "-----------------------------------------------------------------"
echo "Step 1: Setting up Reference Directory"

mkdir -p "$INSAR_DIR"
echo "$LAMBDA" > "$INSAR_DIR/lambda.1.in"

cd "$INSAR_DIR" || exit 1
mkdir -p reference
cd reference || exit 1

echo "Linking Reference SLCs..."
find "$SLC_STACK_PATH/$REF_DATE" -name "$REF_DATE.slc$SLC_SUFFIX*" -exec ln -sf {} . \;

# Handle HDR files
if [ -e "$SLC_STACK_PATH/$REF_DATE/$REF_DATE.slc.hdr" ]; then
    ln -sf "$SLC_STACK_PATH/$REF_DATE/$REF_DATE.slc.hdr" .
elif [ -e "$SLC_STACK_PATH/$REF_DATE/$REF_DATE.hdr" ]; then
    ln -sf "$SLC_STACK_PATH/$REF_DATE/$REF_DATE.hdr" .
fi

# Generate VRT and XML
rm -f "$REF_DATE.slc$SLC_SUFFIX.vrt" "$REF_DATE.slc$SLC_SUFFIX.xml"
gdal_translate -of VRT "$REF_DATE.slc$SLC_SUFFIX" "$REF_DATE.slc$SLC_SUFFIX.vrt"
gdal2isce_xml.py -i "$REF_DATE.slc$SLC_SUFFIX.vrt"

ln -sf "$SLC_STACK_PATH/$REF_DATE/$REF_DATE.slc$SLC_SUFFIX" reference.slc
REF_SLC_PATH="reference/$REF_DATE.slc$SLC_SUFFIX"

# Get Dimensions
cd .. # Back to INSAR_DIR
get_LengthWidth.py "$REF_SLC_PATH"

# --- Step 2: Geometry Setup ---
echo "-----------------------------------------------------------------"
echo "Step 2: Linking Geometry Files"

mkdir -p reference/geom
cd reference/geom || exit 1

GEOM_FILES="los.rdr lat.rdr lon.rdr hgt.rdr"
for GFILE in $GEOM_FILES; do
    SOURCE="$GEOM_PATH/${GFILE}${GEOM_SUFFIX}"
    if ls "${SOURCE}"* 1> /dev/null 2>&1; then
        ln -sf "${SOURCE}"* .
    else
        echo "Error: Geometry file ${SOURCE} not found."
        exit 1
    fi
    # Check for XMLs
    if [ ! -e "${GFILE}${GEOM_SUFFIX}.xml" ]; then
        echo "Error: ${GFILE}${GEOM_SUFFIX}.xml not found."
        exit 1
    fi
done
cd ../.. # Back to INSAR_DIR

# --- Step 3: Integrated Logic (isce2stamps + isce_los2stamps_ENU) ---
echo "-----------------------------------------------------------------"
echo "Step 3: Extracting Geometry & ENU Coefficients"

# Define paths relative to INSAR_DIR (we are currently in INSAR_DIR)
# Source Files
LOS_FILE="reference/geom/los.rdr$GEOM_SUFFIX"
LON_FILE="reference/geom/lon.rdr$GEOM_SUFFIX"
LAT_FILE="reference/geom/lat.rdr$GEOM_SUFFIX"
# Note: In slc_stack mode, the DEM is hgt.rdr
DEM_FILE="reference/geom/hgt.rdr$GEOM_SUFFIX"

# Output Files (Standard StaMPS names)
HEADING_OUT="heading.raw"
INC_OUT="inc_angle.raw"
LON_OUT="lon.raw"
LAT_OUT="lat.raw"
DEM_OUT="dem.raw"

# --- 3.1 Extract DEM (isce2stamps logic) ---
if [ ! -f "$DEM_OUT" ]; then
    echo "  Extracting DEM ($DEM_OUT)..."
    imageMath.py -e="a_0" --a="$DEM_FILE" -o "$DEM_OUT" -s BIL
else
    echo "  $DEM_OUT exists, skipping."
fi

# --- 3.2 Extract Heading (isce2stamps logic) ---
# Logic: -1 * Band_1 - 270 
if [ ! -f "$HEADING_OUT" ]; then
    echo "  Extracting Heading ($HEADING_OUT)..."
    imageMath.py -e='-1*a_1-270' --a="$LOS_FILE" -o "$HEADING_OUT" -s BIL
    
    # Create VRT to handle nodata
    gdal_translate -a_nodata -270 -of VRT "$HEADING_OUT" "$HEADING_OUT.vrt"
    
    # Calculate Mean Heading
    echo "  Calculating mean heading..."
    get_mean_isce.py "$HEADING_OUT" > heading.1.in
else
    echo "  $HEADING_OUT exists, skipping."
fi

# --- 3.3 Extract Incidence Angle (isce2stamps logic) ---
# Logic: Band 0 of LOS file
if [ ! -f "$INC_OUT" ]; then
    echo "  Extracting Incidence Angle ($INC_OUT)..."
    imageMath.py -e="a_0" --a="$LOS_FILE" -o "$INC_OUT" -s BIL
else
    echo "  $INC_OUT exists, skipping."
fi

# --- 3.4 ENU Conversion (INTEGRATED isce_los2stamps_ENU logic) ---
# Previous script called an external tcsh script here. Now we do it inline.
# Formula source: isce_los2stamps_ENU
# East:  sin(rad(a_0))*cos(rad(a_1+90))
# North: sin(rad(a_0))*sin(rad(a_1+90))
# Up:    cos(rad(a_0))
# a_0 is incidence (Band 0), a_1 is azimuth (Band 1)

echo "  Running Integrated ENU Conversion..."

if [ ! -f "e.raw" ]; then
    echo "    Creating East-to-LOS file (e.raw)..."
    imageMath.py --eval='sin(rad(a_0))*cos(rad(a_1+90))' --a="$LOS_FILE" -t FLOAT -s BIL -o "e.raw"
else
    echo "    e.raw exists, skipping."
fi

if [ ! -f "n.raw" ]; then
    echo "    Creating North-to-LOS file (n.raw)..."
    imageMath.py --eval='sin(rad(a_0))*sin(rad(a_1+90))' --a="$LOS_FILE" -t FLOAT -s BIL -o "n.raw"
else
    echo "    n.raw exists, skipping."
fi

if [ ! -f "u.raw" ]; then
    echo "    Creating Up-to-LOS file (u.raw)..."
    imageMath.py --eval='cos(rad(a_0))' --a="$LOS_FILE"  -t FLOAT -s BIL -o "u.raw"
else
    echo "    u.raw exists, skipping."
fi

# --- 3.5 Extract Longitude (isce2stamps logic) ---
if [ ! -f "$LON_OUT" ]; then
    echo "  Extracting Longitude ($LON_OUT)..."
    imageMath.py -e="a_0" --a="$LON_FILE" -o "$LON_OUT" -s BIL
else
    echo "  $LON_OUT exists, skipping."
fi

# --- 3.6 Extract Latitude (isce2stamps logic) ---
if [ ! -f "$LAT_OUT" ]; then
    echo "  Extracting Latitude ($LAT_OUT)..."
    imageMath.py -e="a_0" --a="$LAT_FILE" -o "$LAT_OUT" -s BIL
else
    echo "  $LAT_OUT exists, skipping."
fi

# --- Step 4: Baselines ---
echo "-----------------------------------------------------------------"
echo "Step 4: Processing Baselines"

step_baseline_stack.py -i "$BASELINE_PATH" -m "$REF_DATE"

# Generate SLC list (excluding master)
cd "$SLC_STACK_PATH" || exit 1
ls -d [0-9]*[0-9] | grep -v "$REF_DATE" | awk -F "/" '{print $NF}' > "$WORK_DIR/$INSAR_DIR/slcs.list"

cd "$WORK_DIR/$INSAR_DIR" || exit 1

# --- Step 5: Process Secondary Images & Visualization ---
echo "-----------------------------------------------------------------"
echo "Step 5: Processing Secondaries & Generating Quicklooks"

SLCS=($(cat slcs.list))

for SEC_DATE in "${SLCS[@]}"; do
    echo "Processing: $SEC_DATE"
    
    mkdir -p "$SEC_DATE/coreg_fine"
    cd "$SEC_DATE/coreg_fine" || exit 1
    
    SEC_SLC_SOURCE="$SLC_STACK_PATH/$SEC_DATE/$SEC_DATE.slc$SLC_SUFFIX"
    
    # Links
    ln -sf "${SEC_SLC_SOURCE}"* .
    if [ -e "$SLC_STACK_PATH/$SEC_DATE/$SEC_DATE.slc.hdr" ]; then
        ln -sf "$SLC_STACK_PATH/$SEC_DATE/$SEC_DATE.slc.hdr" .
    elif [ -e "$SLC_STACK_PATH/$SEC_DATE/$SEC_DATE.hdr" ]; then
        ln -sf "$SLC_STACK_PATH/$SEC_DATE/$SEC_DATE.hdr" .
    fi

    # VRT/XML
    rm -f "$SEC_DATE.slc$SLC_SUFFIX.vrt" "$SEC_DATE.slc$SLC_SUFFIX.xml"
    gdal_translate -of VRT "$SEC_DATE.slc$SLC_SUFFIX" "$SEC_DATE.slc$SLC_SUFFIX.vrt"
    gdal2isce_xml.py -i "$SEC_DATE.slc$SLC_SUFFIX.vrt"
    
    ln -sf "$SEC_SLC_SOURCE" coreg.slc
    
    cd ../.. # Back to INSAR_DIR
    
    # Interferogram & Quicklook
    IFG_DIR="$SEC_DATE"
    SAVE_IFG="$IFG_DIR/isce_minrefdem.int"
    REF_SLC_LOC="reference/$REF_DATE.slc$SLC_SUFFIX"
    SEC_SLC_LOC="$SEC_DATE/coreg_fine/$SEC_DATE.slc$SLC_SUFFIX"
    
    if [ "$GENERATE_IFGS" == "y" ]; then
        if [ ! -f "$SAVE_IFG" ]; then
            # Generate Interferogram
            imageMath.py -e='a*conj(b)' --a="$REF_SLC_LOC" --b="$SEC_SLC_LOC" -o "$SAVE_IFG" -t CFLOAT -s BIP >> processing_SM.log
            fixImageXml.py -f -i "$SAVE_IFG"
            
            # Generate Visualization (Quicklook, use mdx.py to visualize)
            LOOKS_IFG="$IFG_DIR/quicklook.int"
            looks.py -i "$SAVE_IFG" -o "$LOOKS_IFG" -r "$RG_LOOKS" -a "$AZ_LOOKS" >> processing_SM.log
            fixImageXml.py -f -i "$LOOKS_IFG"
        fi
    fi
    
    # Final StaMPS Links
    cd "$IFG_DIR" || exit 1
    ln -sf "../reference/reference.slc" reference.slc
    ln -sf "coreg_fine/coreg.slc" secondary.slc
    cd .. # Back to INSAR_DIR
done

echo "-----------------------------------------------------------------"
echo "Finished. Check directory: INSAR_$REF_DATE"

cd "$WORK_DIR" || exit 0