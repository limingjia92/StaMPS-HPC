#!/bin/bash

# ==============================================================================
# Script Name: make_isce_stack_sbas.sh
# Description: Prepares Small Baselines (SBAS) Stack from ISCE SLCs for StaMPS.
#              Replaces 'make_small_baselines_isce'.
#              Process after 'stackSentinel.py'
#
# Original Logic: David Bekaert (2017)
# Refactoring & Optimization: Mingjia Li (2025)
# Language: Bash
# ==============================================================================

# --- Function: Print Usage and Exit ---
print_usage() {
    echo "Usage: $(basename "$0") neighbor_count [slc_stack_path] [stack_geom_path] [slc_stack_baseline_path] [rg_looks] [az_looks]"
    echo ""
    echo "Required Arguments:"
    echo "  neighbor_count          : Number of temporal neighbors to connect (e.g., 2 or 3)"
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

# --- Initialize Variables ---
NEIGHBORS="$1"

# Handle Optional Paths
ARG_SLC="${2:-merged/SLC}"
ARG_GEOM="${3:-merged/geom_reference}"
ARG_BASE="${4:-merged/baselines}"

SLC_STACK_PATH=$(readlink -f "$ARG_SLC")
GEOM_PATH=$(readlink -f "$ARG_GEOM")
BASELINE_PATH=$(readlink -f "$ARG_BASE")

# Visualization Parameters
RG_LOOKS=${5:-20}
AZ_LOOKS=${6:-5}

# --- Hardcoded Parameters ---
SLC_SUFFIX=".full"
GEOM_SUFFIX=".full"
# Sentinel-1 C-band Wavelength (m)
# Note: TerraSAR-X is ~0.031, ALOS-2 is ~0.236. Change this if using other sensors.
LAMBDA=0.0554658 

# Setup Directories
WORK_DIR=$(pwd)
SB_DIR="SMALL_BASELINES"
GENERATE_IFGS="y"

echo "================================================================="
echo "   ISCE to StaMPS SBAS Stack Preparator "
echo "================================================================="
echo "Neighbors   : Connect each date to next $NEIGHBORS dates"
echo "SLC Path    : $SLC_STACK_PATH"
echo "Geom Path   : $GEOM_PATH"
echo "Baselines   : $BASELINE_PATH"
echo "Wavelength  : $LAMBDA"
echo "Quicklooks  : Rg:$RG_LOOKS / Az:$AZ_LOOKS"
echo "================================================================="

# --- Check Inputs ---
if [ ! -d "$SLC_STACK_PATH" ]; then
    echo "Error: SLC path '$SLC_STACK_PATH' does not exist."
    exit 1
fi

if [[ ! "$NEIGHBORS" =~ ^[0-9]+$ ]]; then
    echo "Error: neighbor_count must be an integer."
    exit 1
fi

# --- Get Date List ---
echo "Scanning SLC dates..."
cd "$SLC_STACK_PATH" || exit 1
# Assume folders are named YYYYMMDD. Sort them chronologically.
DATE_LIST=($(ls -d [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9] | sort))
NUM_DATES=${#DATE_LIST[@]}

if [ "$NUM_DATES" -lt 2 ]; then
    echo "Error: Found less than 2 dates. Cannot form baselines."
    exit 1
fi

# Select the FIRST date as the "Super Master" for Geometry Reference
# In SBAS, we need one common grid geometry. We pick the first one.
SUPER_MASTER=${DATE_LIST[0]}
echo "Found $NUM_DATES dates."
echo "Selected Super Master for Geometry: $SUPER_MASTER"

cd "$WORK_DIR" || exit 1
mkdir -p "$SB_DIR"


# --- Step 1: Geometry & Metadata Setup ---
echo "-----------------------------------------------------------------"
echo "Step 1: Setting up Geometry (Super Master: $SUPER_MASTER)"


echo "$LAMBDA" > "$SB_DIR/lambda.1.in"
cd "$SB_DIR" || exit 1

# 1.1 Link Geometry Files
mkdir -p geometry_ref
cd geometry_ref || exit 1

echo "Linking Geometry files..."
GEOM_FILES="los.rdr lat.rdr lon.rdr hgt.rdr"
for GFILE in $GEOM_FILES; do
    SOURCE="$GEOM_PATH/${GFILE}${GEOM_SUFFIX}"
    if ls "${SOURCE}"* 1> /dev/null 2>&1; then
        ln -sf "${SOURCE}"* .
    else
        echo "Error: Geometry file ${SOURCE} not found."
        exit 1
    fi
    # Check XML
    if [ ! -e "${GFILE}${GEOM_SUFFIX}.xml" ]; then
        echo "Error: ${GFILE}${GEOM_SUFFIX}.xml not found."
        exit 1
    fi
done

cd .. # Back to SMALL_BASELINES root

# --- Step 2: Integrated Logic (isce2stamps + isce_los2stamps_ENU) ---
echo "-----------------------------------------------------------------"
echo "Step 2: Extracting Geometry & ENU Coefficients"

# Source paths
LOS_FILE="geometry_ref/los.rdr$GEOM_SUFFIX"
LON_FILE="geometry_ref/lon.rdr$GEOM_SUFFIX"
LAT_FILE="geometry_ref/lat.rdr$GEOM_SUFFIX"
# In slc_stack mode, DEM is hgt.rdr
DEM_FILE="geometry_ref/hgt.rdr$GEOM_SUFFIX"

# Output files
HEADING_OUT="heading.raw"
INC_OUT="inc_angle.raw"
LON_OUT="lon.raw"
LAT_OUT="lat.raw"
DEM_OUT="dem.raw"

# Extract DEM
if [ ! -f "$DEM_OUT" ]; then
    imageMath.py -e="a_0" --a="$DEM_FILE" -o "$DEM_OUT" -s BIL
fi

# Extract Heading
if [ ! -f "$HEADING_OUT" ]; then
    imageMath.py -e='-1*a_1-270' --a="$LOS_FILE" -o "$HEADING_OUT" -s BIL
    gdal_translate -a_nodata -270 -of VRT "$HEADING_OUT" "$HEADING_OUT.vrt"
    get_mean_isce.py "$HEADING_OUT" > heading.1.in
fi

# Extract Incidence Angle
if [ ! -f "$INC_OUT" ]; then
    imageMath.py -e="a_0" --a="$LOS_FILE" -o "$INC_OUT" -s BIL
fi

# Extract Lon/Lat
if [ ! -f "$LON_OUT" ]; then
    imageMath.py -e="a_0" --a="$LON_FILE" -o "$LON_OUT" -s BIL
fi
if [ ! -f "$LAT_OUT" ]; then
    imageMath.py -e="a_0" --a="$LAT_FILE" -o "$LAT_OUT" -s BIL
fi

# 1.3 ENU Conversion (Internal isce_los2stamps_ENU Logic)
echo "Calculating ENU vectors..."
if [ ! -f "e.raw" ]; then
    imageMath.py --eval='sin(rad(a_0))*cos(rad(a_1+90))' --a="$LOS_FILE" -t FLOAT -s BIL -o "e.raw"
fi
if [ ! -f "n.raw" ]; then
    imageMath.py --eval='sin(rad(a_0))*sin(rad(a_1+90))' --a="$LOS_FILE" -t FLOAT -s BIL -o "n.raw"
fi
if [ ! -f "u.raw" ]; then
    imageMath.py --eval='cos(rad(a_0))' --a="$LOS_FILE"  -t FLOAT -s BIL -o "u.raw"
fi

# 1.4 Get Dimensions
# We need to link the Super Master SLC temporarily to get width/length
mkdir -p tmp_ref
cd tmp_ref
ln -sf "$SLC_STACK_PATH/$SUPER_MASTER/$SUPER_MASTER.slc$SLC_SUFFIX" .
ln -sf "$SLC_STACK_PATH/$SUPER_MASTER/$SUPER_MASTER.slc.hdr" . 2>/dev/null
ln -sf "$SLC_STACK_PATH/$SUPER_MASTER/$SUPER_MASTER.hdr" . 2>/dev/null

# Generate XML if needed
rm -f *.vrt *.xml
gdal_translate -of VRT "$SUPER_MASTER.slc$SLC_SUFFIX" "$SUPER_MASTER.slc$SLC_SUFFIX.vrt"
gdal2isce_xml.py -i "$SUPER_MASTER.slc$SLC_SUFFIX.vrt"

cd .. # Back to SMALL_BASELINES
get_LengthWidth.py "tmp_ref/$SUPER_MASTER.slc$SLC_SUFFIX"
rm -rf tmp_ref # Cleanup

# --- Step 3: Baselines ---
echo "-----------------------------------------------------------------"
echo "Step 3: Processing Baselines"
# We calculate baselines relative to our Super Master ($SUPER_MASTER)

step_baseline_stack.py -i "$BASELINE_PATH" -m "$SUPER_MASTER"

# --- Step 4: Interferogram Generation (Auto-Pairing Loop) ---
echo "-----------------------------------------------------------------"
echo "Step 4: Generating Interferograms (Neighbors: $NEIGHBORS)"

# Re-read date list into array (Bash indices start at 0)
DATES=("${DATE_LIST[@]}")
COUNT=${#DATES[@]}

# Loop through dates
for (( i=0; i<$COUNT; i++ )); do
    MASTER_DATE=${DATES[$i]}
    
    # Inner loop for neighbors
    for (( j=1; j<=$NEIGHBORS; j++ )); do
        TARGET_IDX=$((i + j))
        
        # Check bound
        if [ $TARGET_IDX -lt $COUNT ]; then
            SLAVE_DATE=${DATES[$TARGET_IDX]}
            
            # === Process Pair: MASTER_DATE - SLAVE_DATE ===
            PAIR_DIR="${MASTER_DATE}_${SLAVE_DATE}"
            echo "Processing Pair: $PAIR_DIR"
            
            mkdir -p "$PAIR_DIR"
            cd "$PAIR_DIR" || exit 1
            
            # 3.1 Define Source Paths
            M_SLC_SRC="$SLC_STACK_PATH/$MASTER_DATE/$MASTER_DATE.slc$SLC_SUFFIX"
            S_SLC_SRC="$SLC_STACK_PATH/$SLAVE_DATE/$SLAVE_DATE.slc$SLC_SUFFIX"
            
            # 3.2 Link SLCs
            # Link Master
            ln -sf "$M_SLC_SRC" reference.slc
            # Link Slave
            ln -sf "$S_SLC_SRC" secondary.slc
            
            # Ensure XML/VRT exists for ISCE (imageMath needs them)
            # Helper function to ensure XML exists
            ensure_xml() {
                local DATE=$1
                local SLC_PATH=$2
                local LOCAL_LINK=$3
                
                # If source XML doesn't exist, we must generate it locally
                # Check source .xml
                if [ ! -f "${SLC_PATH}.xml" ]; then
                    # We need to fake it locally
                    # This part is tricky if source folder is read-only
                    # Assuming we can make VRT locally based on the link
                    gdal_translate -of VRT "$LOCAL_LINK" "$LOCAL_LINK.vrt"
                    gdal2isce_xml.py -i "$LOCAL_LINK.vrt"
                else
                   # Link the XML from source if it exists
                   ln -sf "${SLC_PATH}.xml" "${LOCAL_LINK}.xml"
                   ln -sf "${SLC_PATH}.vrt" "${LOCAL_LINK}.vrt"
                fi
            }
            
            ensure_xml "$MASTER_DATE" "$M_SLC_SRC" "reference.slc"
            ensure_xml "$SLAVE_DATE"  "$S_SLC_SRC" "secondary.slc"
            
            # 3.3 Compute Interferogram
            SAVE_IFG="isce_minrefdem.int"
            
            if [ "$GENERATE_IFGS" == "y" ]; then
                if [ ! -f "$SAVE_IFG" ]; then
                    # imageMath
                    imageMath.py -e='a*conj(b)' --a="reference.slc" --b="secondary.slc" -o "$SAVE_IFG" -t CFLOAT -s BIP >> ../processing_SB.log
                    fixImageXml.py -f -i "$SAVE_IFG"
                    
                    # 3.4 Visualization (Quicklook)
                    echo "  Generating Quicklook..."
                    LOOKS_IFG="quicklook.int"
                    looks.py -i "$SAVE_IFG" -o "$LOOKS_IFG" -r "$RG_LOOKS" -a "$AZ_LOOKS" >> ../processing_SB.log
                    fixImageXml.py -f -i "$LOOKS_IFG"
                fi
            fi
            # --- Cleanup: Remove temporary SLC links ---
            rm -f reference.slc secondary.slc
            rm -f reference.slc.xml reference.slc.vrt secondary.slc.xml secondary.slc.vrt
            
            cd .. # Back to SMALL_BASELINES
        fi
    done
done

# --- Step 5: SLC Directory Setup ---
echo "-----------------------------------------------------------------"
echo "Step 5: Setting up individual SLC directories"
# This creates a folder for each date and links the SLC, 
# mimicking the structure expected by some StaMPS steps.

for DATE in "${DATES[@]}"; do
    if [ ! -d "$DATE" ]; then
        mkdir "$DATE"
    fi
    
    cd "$DATE" || exit 1
    
    # Define source SLC path
    SRC_SLC="$SLC_STACK_PATH/$DATE/$DATE.slc$SLC_SUFFIX"
    
    # Create the link: Date.slc -> /path/to/Date.slc.full
    # Using relative link name as requested: "20241128.slc"
    TARGET_NAME="${DATE}.slc"
    
    # Check if link already exists
    if [ ! -e "$TARGET_NAME" ]; then
        echo "Linking SLC for $DATE..."
        ln -sf "$SRC_SLC" "$TARGET_NAME"
        
        # Optional: Link XML/VRT if they exist in source (Helpful for StaMPS)
        if [ -f "${SRC_SLC}.xml" ]; then ln -sf "${SRC_SLC}.xml" "${TARGET_NAME}.xml"; fi
        if [ -f "${SRC_SLC}.vrt" ]; then ln -sf "${SRC_SLC}.vrt" "${TARGET_NAME}.vrt"; fi
        if [ -f "${SRC_SLC}.hdr" ]; then ln -sf "${SRC_SLC}.hdr" "${TARGET_NAME}.hdr"; fi
    fi
    
    cd .. # Back to SMALL_BASELINES
done

echo "-----------------------------------------------------------------"
echo "Finished. Output directory: $WORK_DIR/$SB_DIR"
echo "Total dates: $COUNT"
echo "Strategy: $NEIGHBORS sequential connections."

echo "-----------------------------------------------------------------"
echo "Finished. Output directory: $WORK_DIR/$SB_DIR"
echo "Total dates: $COUNT"
echo "Strategy: $NEIGHBORS sequential connections."

cd "$WORK_DIR" || exit 0