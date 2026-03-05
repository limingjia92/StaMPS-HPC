#!/bin/bash

# ==============================================================================
# Script Name: prep_stamps_gamma.sh
# Description: Prepares GAMMA output for StaMPS processing with parallel patching.
#
# Original Author: Andy Hooper, December 2012
# Optimization & Parallelization: Mingjia Li, November 2025
# 
# Note on Parallelization:
# This script supports parallel processing via 'jobmaxnum'. However, if using
# mechanical hard drives (HDD), increasing job concurrency may NOT improve 
# speed (and could be slower) due to I/O read/write bottlenecks.
# ==============================================================================

# --- Function: Print Usage and Exit ---
print_usage() {
    echo "Usage: $(basename "$0") reference_date data_dir work_dir [da_thresh] [rg_patches] [az_patches] [jobmaxnum] [rg_overlap] [az_overlap] [maskfile]"
    echo ""
    echo "Arguments:"
    echo "  reference_date    : Reference image date (e.g., 20191126)"
    echo "  data_dir          : Path to input data directory (must follow structure below)"
    echo "  work_dir          : Path to output working directory"
    echo "  da_thresh         : Amplitude dispersion threshold (lower = stricter & fewer points)"
    echo "                      (Default: 0.4 for PS, 0.6 for SB)"
    echo "  rg_patches        : Number of patches in Range direction (Default: 1)"
    echo "  az_patches        : Number of patches in Azimuth direction (Default: 1)"
    echo "  jobmaxnum         : Max parallel jobs (Default: 1)"
    echo "                      *WARNING: High concurrency on HDDs may cause I/O bottlenecks."
    echo "  rg_overlap        : Overlapping pixels in Range (Default: 400)"
    echo "  az_overlap        : Overlapping pixels in Azimuth (Default: 400)"
    echo "  maskfile          : (Optional) Binary mask file (0=include, 1=exclude)"
    echo ""
    echo "Selective Patch Processing (Advanced):"
    echo "  If you wish to process only a specific subset of patches:"
    echo "  1. Use the MATLAB tool 'plot_patch_kml_gamma' to visualize patch layout and identify IDs."
    echo "  2. Create or edit a file named 'patch.list' inside 'work_dir'."
    echo "  3. List only the Patch directories you want to process (e.g., PATCH_1, PATCH_5)."
    echo "  4. Run this script. It will detect the existing 'patch.list' and ONLY process the listed patches."
    echo ""
    echo "Expected Data Directory Structure:"
    echo "  (Note: For SB mode, 'rslc' and 'diff0' must be inside 'SMALL_BASELINES/')"
    echo ""
    echo "    [SMALL_BASELINES/]rslc/*.rslc        # SLC files"
    echo "    [SMALL_BASELINES/]rslc/*.rslc.par    # SLC parameter files"
    echo "    [SMALL_BASELINES/]diff0/*.diff       # Differential interferograms"
    echo "    [SMALL_BASELINES/]diff0/*.base       # Baseline files"
    echo "    geo/*diff_par                        # Diff parameter file"
    echo "    geo/*dem.rdc                         # Radar coordinates DEM"
    echo "    geo/YYYYMMDD.lon                     # Longitude (Radar Coords)"
    echo "    geo/YYYYMMDD.lat                     # Latitude (Radar Coords)"
    echo ""
    exit 1
}

# --- Check Arguments ---
if [ $# -lt 3 ]; then
    print_usage
fi

# --- Initialize Variables ---
REFERENCE_DATE=$1
DATA_DIR=$2
WORK_DIR=$3

# Optional Arguments with Defaults (handling '-' as default placeholder)
DA_THRESH=${4:--}
PATCHES_RG=${5:--}; [[ "$PATCHES_RG" == "-" ]] && PATCHES_RG=1
PATCHES_AZ=${6:--}; [[ "$PATCHES_AZ" == "-" ]] && PATCHES_AZ=1
MAX_JOBS=${7:--};   [[ "$MAX_JOBS" == "-" ]] && MAX_JOBS=1
OVERLAP_RG=${8:--}; [[ "$OVERLAP_RG" == "-" ]] && OVERLAP_RG=400
OVERLAP_AZ=${9:--}; [[ "$OVERLAP_AZ" == "-" ]] && OVERLAP_AZ=400
MASK_FILE=${10:-""}

TIME_START=$(date +"%s")

echo "================================================================="
echo "   StaMPS Pre-processor (GAMMA Interface) - HPC Optimized"
echo "   Base Logic: Andy Hooper (2012)"
echo "   Parallelization: Mingjia Li (2025)"
echo "================================================================="
echo "Reference Date: $REFERENCE_DATE"
echo "Data Dir   : $DATA_DIR"
echo "Work Dir   : $WORK_DIR"
echo "Patches    : $PATCHES_RG (Range) x $PATCHES_AZ (Azimuth)"
echo "Max Jobs   : $MAX_JOBS"

# Check Mask File
if [ -n "$MASK_FILE" ]; then
    if [ ! -e "$MASK_FILE" ]; then
        echo "Error: Mask file $MASK_FILE does not exist."
        exit 2
    fi
    echo "Mask File  : $MASK_FILE"
fi

# Ensure Work Directory
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit 1

# --- Detect Mode (PS or SB) ---
if [ -d "$DATA_DIR/SMALL_BASELINES" ]; then
    echo "Mode       : Small Baselines (SB)"
    SB_FLAG=1
    # Find the reference parameter file in SB structure
    RSC_FILE=$(ls "$DATA_DIR"/SMALL_BASELINES/*/"$REFERENCE_DATE".*slc.par | head -n 1)
else
    echo "Mode       : Persistent Scatterers (PS)"
    SB_FLAG=0
    # Find the reference parameter file in PS structure
    RSC_FILE=$(ls "$DATA_DIR"/*slc/"$REFERENCE_DATE".*slc.par | head -n 1)
fi

echo "Param File : $RSC_FILE"

# Set Threshold if not provided
if [ "$DA_THRESH" == "-" ]; then
    if [ $SB_FLAG -eq 1 ]; then
        DA_THRESH=0.6
    else
        DA_THRESH=0.4
    fi
fi
echo "Amp Thresh : $DA_THRESH"

# --- Extract Image Metadata ---
IMG_LENGTH=$(gawk '/azimuth_lines/ {print $2}' < "$RSC_FILE")
IMG_WIDTH=$(gawk '/range_samples/ {print $2}' < "$RSC_FILE")
IMG_FORMAT=$(gawk '/image_format/ {print $2}' < "$RSC_FILE")

echo "Dimensions : $IMG_WIDTH (Width) x $IMG_LENGTH (Length)"

if [ "$IMG_FORMAT" == "FCOMPLEX" ]; then
    PRECISION="f"
elif [ "$IMG_FORMAT" == "SCOMPLEX" ]; then
    PRECISION="s"
else
    echo "Error: Image format is $IMG_FORMAT. Must be FCOMPLEX or SCOMPLEX."
    exit 1
fi
echo "Precision  : $PRECISION"

# --- Initial Setup ---
echo "gamma" > processor.txt

# Run Matlab parameter initialization
echo "Running Matlab parameter initialization..."
if [ $SB_FLAG -eq 1 ]; then
    matlab -nojvm -nosplash -nodisplay -r "sb_parms_initial; exit" > sb_parms_initial.log
else
    matlab -nojvm -nosplash -nodisplay -r "ps_parms_initial; exit" > ps_parms_initial.log
fi

# Save metadata for later steps
echo "$IMG_WIDTH" > "$WORK_DIR/width.txt"
echo "$IMG_LENGTH" > "$WORK_DIR/len.txt"
echo "$RSC_FILE" > "$WORK_DIR/rsc.txt"

# --- Calibrate Amplitudes (calamp) ---
echo "-----------------------------------------------------------------"
echo "Step 1: Amplitude Calibration"

rm -f "$WORK_DIR/calamp.in" 2>/dev/null

if [ $SB_FLAG -eq 1 ]; then
    ls "$DATA_DIR"/SMALL_BASELINES/*/*.*slc >> "$WORK_DIR/calamp.in"
    SEL_FILE="$WORK_DIR/selsbc.in"
else
    ls "$DATA_DIR"/*slc/*.*slc >> "$WORK_DIR/calamp.in"
    SEL_FILE="$WORK_DIR/selpsc.in"
fi

if [ -s "$WORK_DIR/calamp.out" ]; then
    echo "Skipping calamp (calamp.out exists)."
else
    # Command syntax: calamp <input_list> <width> <output> <type> <1=swap byte> <maskfile>
    echo "Running calamp..."
    
    # Optimization: Set OpenMP threads to 4 to speed up reading SLCs
    export OMP_NUM_THREADS=4 
    calamp calamp.in "$IMG_WIDTH" calamp.out "$PRECISION" 1 "$MASK_FILE" > calamp.log
    sort calamp.out -o calamp.out
fi

# Prepare selection input file
rm -f "$SEL_FILE" 2>/dev/null
echo "$DA_THRESH" > "$SEL_FILE"
echo "$IMG_WIDTH" >> "$SEL_FILE"
cat "$WORK_DIR/calamp.out" >> "$SEL_FILE"

# --- Patch Division ---
echo "-----------------------------------------------------------------"
echo "Step 2: Dividing Image into Patches"

WIDTH_P=$((IMG_WIDTH / PATCHES_RG))
LENGTH_P=$((IMG_LENGTH / PATCHES_AZ))
echo "Patch Size : $WIDTH_P x $LENGTH_P pixels (approx)"

# Check if patch.list exists to determine mode
if [ -f "patch.list" ]; then
    echo "Found existing patch.list. Mode: Check and Fill Missing Directories."
    EXISTING_LIST_MODE=1
else
    echo "No patch.list found. Mode: Create New List and Directories."
    rm -rf patch.list
    EXISTING_LIST_MODE=0
fi

PATCH_COUNT=0

for ((irg=1; irg<=PATCHES_RG; irg++)); do
    for ((iaz=1; iaz<=PATCHES_AZ; iaz++)); do
        PATCH_COUNT=$((PATCH_COUNT + 1))
        CURRENT_PATCH_DIR="PATCH_$PATCH_COUNT"
        
        # --- Decision Logic: Should we process this patch? ---
        PROCESS_THIS_PATCH=0
        
        if [ $EXISTING_LIST_MODE -eq 1 ]; then
            # 1. Existing List Mode
            # Check if this patch name is in the existing file (exact match)
            if grep -Fxq "$CURRENT_PATCH_DIR" patch.list; then
                if [ -d "$CURRENT_PATCH_DIR" ]; then
                    # Dir exists, skip
                    # echo "  $CURRENT_PATCH_DIR exists. Skipping."
                    continue
                else
                    # In list but dir missing, need to create
                    echo "  $CURRENT_PATCH_DIR found in list but missing on disk. Generating..."
                    PROCESS_THIS_PATCH=1
                fi
            else
                # Not in the existing list (user might have removed it), skip
                continue
            fi
        else
            # 2. New List Mode (Original Behavior)
            echo "$CURRENT_PATCH_DIR" >> patch.list
            PROCESS_THIS_PATCH=1
            # Check dir to avoid error, though usually we'd overwrite
            if [ -d "$CURRENT_PATCH_DIR" ]; then
                # Optional: warn if overwriting, currently we just regenerate contents
                : 
            fi
        fi

        # --- Generation Logic (Only runs if PROCESS_THIS_PATCH=1) ---
        if [ $PROCESS_THIS_PATCH -eq 1 ]; then
            
            # Calculate Range Indices
            START_RG1=$((WIDTH_P * (irg - 1) + 1))
            START_RG=$((START_RG1 - OVERLAP_RG))
            [[ $START_RG -lt 1 ]] && START_RG=1
            
            END_RG1=$((WIDTH_P * irg))
            END_RG=$((END_RG1 + OVERLAP_RG))
            [[ $END_RG -gt $IMG_WIDTH ]] && END_RG=$IMG_WIDTH
            
            # Calculate Azimuth Indices
            START_AZ1=$((LENGTH_P * (iaz - 1) + 1))
            START_AZ=$((START_AZ1 - OVERLAP_AZ))
            [[ $START_AZ -lt 1 ]] && START_AZ=1
            
            END_AZ1=$((LENGTH_P * iaz))
            END_AZ=$((END_AZ1 + OVERLAP_AZ))
            [[ $END_AZ -gt $IMG_LENGTH ]] && END_AZ=$IMG_LENGTH

            # Create Patch Directory if missing
            if [ ! -d "$CURRENT_PATCH_DIR" ]; then
                mkdir "$CURRENT_PATCH_DIR"
            fi

            # Write patch parameters
            cd "$CURRENT_PATCH_DIR" || exit 1
            
            # patch.in (with overlaps)
            printf "%s\n%s\n%s\n%s\n" "$START_RG" "$END_RG" "$START_AZ" "$END_AZ" > patch.in
            
            # patch_noover.in (without overlaps)
            printf "%s\n%s\n%s\n%s\n" "$START_RG1" "$END_RG1" "$START_AZ1" "$END_AZ1" > patch_noover.in
            
            cd "$WORK_DIR" || exit 1
        fi
    done
done

# --- Prepare Global Input Files ---
echo "-----------------------------------------------------------------"
echo "Step 3: Preparing Global Input Files"

# 1. pscphase.in
echo "$IMG_WIDTH" > pscphase.in
if [ $SB_FLAG -eq 1 ]; then
    ls "$DATA_DIR"/SMALL_BASELINES/*/*.diff >> pscphase.in
else
    ls "$DATA_DIR"/diff0/*.diff >> pscphase.in
fi

# 2. pscdem.in
echo "$IMG_WIDTH" > pscdem.in
ls "$DATA_DIR"/geo/*dem.rdc >> pscdem.in

# 3. psclonlat.in (Mandatory check)
LON_FILE=$(ls -1 "$DATA_DIR"/geo/*.lon 2>/dev/null | head -1)
LAT_FILE=$(ls -1 "$DATA_DIR"/geo/*.lat 2>/dev/null | head -1)

if [[ -n "$LON_FILE" && -n "$LAT_FILE" ]]; then
    echo "$IMG_WIDTH" > psclonlat.in
    echo "$LON_FILE" >> psclonlat.in
    echo "$LAT_FILE" >> psclonlat.in
    echo "Geo Check  : .lon and .lat files found. Using psclonlat."
else
    echo "ERROR: .lon or .lat files not found in $DATA_DIR/geo/"
    echo "       This script enforces 'psclonlat' and does not support pt2geo fallback."
    exit 1
fi

# 4. pscheading.in
ls "$DATA_DIR"/*/*dem_par > pscheading.in 2>/dev/null
ls "$DATA_DIR"/*/lv_theta >> pscheading.in 2>/dev/null
ls "$DATA_DIR"/*/lv_phi >> pscheading.in 2>/dev/null

# --- Parallel Execution ---
echo "-----------------------------------------------------------------"
echo "Step 4: Parallel Processing (Max Jobs: $MAX_JOBS)"

BYTESWAP=1
LIST_FILE="patch.list"
TOTAL_PATCHES=$(wc -l < "$WORK_DIR/$LIST_FILE")

# FIFO Setup for concurrency control
FIFO_NAME="tempfifo_$$"
trap "rm -f $FIFO_NAME; exit 1" INT TERM
mkfifo "$FIFO_NAME"
exec 1001<>"$FIFO_NAME"
rm -f "$FIFO_NAME"

# Pre-fill FIFO tokens
for ((i=0; i<MAX_JOBS; i++)); do echo >&1001; done

echo "Processing $TOTAL_PATCHES patches..."

for ((i=1; i<=TOTAL_PATCHES; ++i)); do
    PATCH_NAME=$(sed -n "${i}p" "$WORK_DIR/$LIST_FILE")
    
    # Wait for a token
    read -u1001
    
    {
        cd "$WORK_DIR/$PATCH_NAME" || exit 1
        LOG_FILE="mt_${PATCH_NAME}.log"
        rm -f "$LOG_FILE"
        
        echo "Processing $PATCH_NAME..." | tee -a "$LOG_FILE"
        
        # 1. Select Candidates (selpsc/selsbc)
        if [ $SB_FLAG -eq 1 ]; then
             CMD="selsbc_patch $WORK_DIR/selsbc.in patch.in pscands.1.ij pscands.1.da mean_amp.flt $PRECISION $BYTESWAP $MASK_FILE"
        else
             CMD="selpsc_patch $WORK_DIR/selpsc.in patch.in pscands.1.ij pscands.1.da mean_amp.flt $PRECISION $BYTESWAP $MASK_FILE"
        fi
        
        # Remove trailing space if MASK_FILE is empty
        CMD=$(echo "$CMD" | xargs)
        
        echo "$CMD" >> "$LOG_FILE"
        $CMD >> "$LOG_FILE" 2>&1
        
        # Calculate PS number
        if [ -f "pscands.1.ij" ]; then
            ps_num=$(wc -l < pscands.1.ij)
        else
            ps_num=0
        fi
        
        echo "" >> "$LOG_FILE"
        echo "Total $ps_num PS candidates extracted" >> "$LOG_FILE"

        # 2. Extract Lon/Lat (psclonlat) - Mandatory
        echo "" >> "$LOG_FILE"
        echo "psclonlat $WORK_DIR/psclonlat.in pscands.1.ij pscands.1.ll" >> "$LOG_FILE"
        psclonlat "$WORK_DIR/psclonlat.in" pscands.1.ij pscands.1.ll >> "$LOG_FILE" 2>&1

        # 3. Extract Height (pscdem)
        echo "" >> "$LOG_FILE"
        echo "pscdem $WORK_DIR/pscdem.in pscands.1.ij pscands.1.hgt" >> "$LOG_FILE"
        pscdem "$WORK_DIR/pscdem.in" pscands.1.ij pscands.1.hgt >> "$LOG_FILE" 2>&1

        # 4. Extract Phase (pscphase)
        echo "" >> "$LOG_FILE"
        echo "pscphase $WORK_DIR/pscphase.in pscands.1.ij pscands.1.ph" >> "$LOG_FILE"
        pscphase "$WORK_DIR/pscphase.in" pscands.1.ij pscands.1.ph >> "$LOG_FILE" 2>&1

        # 5. Extract Heading and Incidence angles (pscheading) - Optional
        line_count=$(wc -l < "$WORK_DIR/pscheading.in")
        if [ "$line_count" -eq 3 ]; then
            echo "" >> "$LOG_FILE"
            echo "pscheading $WORK_DIR/pscheading.in pscands.1.ll pscands.1.head pscands.1.inc" >> "$LOG_FILE"
            pscheading "$WORK_DIR/pscheading.in" pscands.1.ll pscands.1.head pscands.1.inc >> "$LOG_FILE" 2>&1
        else
            echo "" >> "$LOG_FILE"
            echo "Skipping pscheading, need result file from look_vector" >> "$LOG_FILE"
        fi
        
        # Return token
        echo >&1001
        echo "$PATCH_NAME complete."
    } &
done

wait

# Cleanup FIFO
exec 1001>&-

# --- Final Reporting ---
TIME_END=$(date +"%s")
DURATION=$((TIME_END - TIME_START))

# Calculate Hours, Minutes
HOURS=$((DURATION / 3600))
MINUTES=$(( (DURATION % 3600) / 60 ))

echo "-----------------------------------------------------------------"
echo "Processing complete."

# Format: Xh Ym 
echo "Total Time Used: ${HOURS}h ${MINUTES}m "

exit 0