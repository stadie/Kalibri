#!/bin/bash

# $Id: createJECValidationHtmlPage.sh,v 1.2 2010/04/23 08:56:39 mschrode Exp $

#  This script creates an html webpage listing JEC validation
#  plots.
#
#  It takes the following parameters as input:
#
#   1: MODE
#   Indicates whether validation plots from Kalibri or JetMET
#   calibration are to be shown. Possible values are 'kalibri'
#   and 'jetmet'.
#
#   2: PLOTS_PATH
#   Path to an .tar archive file storing the validation plots.
#   The files are expected to be in eps or pdf format and to conform
#   to the Kalibri naming conventions (compare the global
#   variables ID_* below):
#    - File names of plots showing quantities vs genJetPt must
#      contain 'ResponseVsGenJetPt'.
#    - File names of plots showing quantities vs eta must
#      contain 'ResponseVsEta'.
#    - File names of plots showing the mean response must
#      contain 'GaussFitMean'.
#    - File names of plots showing the mean resolution must
#      contain 'GaussFitWidth'.
#    - File names of plots having a zoomed y axis must
#      *end with* 'zoom.eps' (or 'zoom.pdf').
#   All the above should be met automatically by the Kalibri
#   ControlPlots class.
#
#   3: ID
#   String defining the validated JEC e.g. "Summer09 7 TeV
#   (rereco)". The id will be used to generate the title of
#   the webpage as well as the name of the directory storing
#   the plots and html file. For the directory name, all
#   characters that are not letters or numbers are replaced
#   by an underscore i.e. space, brackets, ... can be used
#   to have a nice webpage title! Note: the jet type needs
#   not to be specified; it is the next argument.
#
#   4: JEC_TAG
#   Official tag of the JEC specified in the correction module i.e.
#   'JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_ReReco332_cff'
#
#   5: DATASET
#   Name of the dataset used to produce the validation plots i.e.
#   '/QCDFlat_Pt15to3000/Summer09-MC_31X_V9_7TeV-v1/GEN-SIM-RECO'
#
#   6: JET_ALGO
#   Used jet algorithm. Possible algorithms are
#    - 'ak5' for anti-kt 5
#    - 'sc5' for sis-cone 5
#
#   7: JET_TYPE
#   Used jet type. Possible types are
#    - 'calo' for calorimeter jets
#    - 'jpt' for jet-plus-track jets
#    - 'pf' for particle-flow jets
#
#   Author: Matthias Schroeder
#   Date: Tue Apr 20 21:17:49 CEST 2010




# ----- Global variables --------------------------------------

# Name of produced html fil
FILE_OUT="index.html"

# Naming conventions for plots
ID_GenPt="GenJetPt"
ID_Eta="Eta"
ID_RVs="ResponseVs"
ID_RVsGenPt="${ID_RVs}${ID_GenPt}"
ID_RVsEta="${ID_RVs}${ID_Eta}"
ID_Response="GaussFitMean"
ID_Resolution="GaussFitWidth"
ID_ZOOM="zoom"

# Titles for different categories of validation plots
TITLE_ResponseVsGenPt="Response vs genJetPt in bins of eta"
TITLE_ResponseVsEta="Response vs eta in bins of genJetPt"
TITLE_ResolutionVsGenPt="Resolution vs genJetPt in bins of eta"
TITLE_ResolutionVsEta="Resolution vs eta in bins of genJetPt"




# ----- Read parameters ---------------------------------------

# Kalibri or JetMET JEC validation
MODE=($1)
CORRECT_INPUT=0
until [[ ${CORRECT_INPUT} -eq 1 ]]; do
    if [[ ${MODE} == "kalibri" ]]; then
	CORRECT_INPUT=1
	JET_TYPE="Kalibri"
    elif [[ ${MODE} == "jetmet" ]]; then
	CORRECT_INPUT=1
	JET_TYPE="JetMET"
    else	
	echo -n "Enter mode ('kalibri','jetmet'): "
	read MODE
    fi
done    


# Name of archive file storing validation plots
PLOTS_PATH=($2)
CORRECT_INPUT=0
until [[ ${CORRECT_INPUT} -eq 1 ]]; do
    if [[ -e ${PLOTS_PATH} ]]; then
	if [[ ${PLOTS_PATH} == *.tar ]]; then
	    CORRECT_INPUT=1
	fi
    fi
    if [[ ${CORRECT_INPUT} -eq 0 ]]; then
	echo -n "Enter name of tar-archive file storing validation plots: "
	read PLOTS_PATH
    fi
done    


# Id i.e. some unofficial name of the JEC
ID=($3)
while [[ -z ${ID} ]]; do
    echo -n "Enter id of JEC (e.g. Summer09 7TeV): "
    read ID
done


# JEC tag name
JEC_TAG=($4)
while [[ -z ${JEC_TAG} ]]; do
    echo -n "Enter JEC tag name: "
    read JEC_TAG
done


# Dataset name
DATASET=($5)
while [[ -z ${DATASET} ]]; do
    echo -n "Enter dataset name: "
    read DATASET
done


# Jet algorithm ('ak5','sc5')
JET_ALGO=($6)
CORRECT_INPUT=0
until [[ ${CORRECT_INPUT} -eq 1 ]]; do
    if [[ ${JET_ALGO} == "ak5" ]]; then
	CORRECT_INPUT=1
	JET_ALGO="AK5"
    elif [[ ${JET_ALGO} == "sc5" ]]; then
	CORRECT_INPUT=1
	JET_ALGO="SC5"
    else	
	echo -n "Enter jet algorithm ('ak5','sc5'): "
	read JET_ALGO
    fi
done    


# Jet type ('calo','jpt','pf')
JET_TYPE=($7)
CORRECT_INPUT=0
until [[ ${CORRECT_INPUT} -eq 1 ]]; do
    if [[ ${JET_TYPE} == "calo" ]]; then
	CORRECT_INPUT=1
	JET_TYPE="Calo"
    elif [[ ${JET_TYPE} == "jpt" ]]; then
	CORRECT_INPUT=1
	JET_TYPE="JPT"
    elif [[ ${JET_TYPE} == "pf" ]]; then
	CORRECT_INPUT=1
	JET_TYPE="PF"
    else	
	echo -n "Enter jet type ('calo','jpt','pf'): "
	read JET_TYPE
    fi
done    




echo "Creating validation webpage for '${ID}'"
echo ""
echo " Mode               :  $MODE"
echo " Archive file name  :  $PLOTS_PATH"
echo " Id                 :  $ID"
echo " JEC tag name       :  $JEC_TAG"
echo " Dataset name       :  $DATASET"
echo " Jet algorithm      :  $JET_ALGO"
echo " Jet type           :  $JET_TYPE"
echo ""




# ----- Prepare plots -----------------------------------------

# Setup html directory i.e. removing non-characters and
# non-numbers from id
DIR_OUT=(${ID//[^a-zA-Z0-9]/_})
DIR_OUT=(${DIR_OUT//_[_*]/_})
DIR_OUT=(${DIR_OUT%_})
DIR_OUT=(${DIR_OUT}_${JET_ALGO})
DIR_OUT=(${DIR_OUT}_${JET_TYPE})
echo "Preparing working directory '$DIR_OUT'"
mkdir $DIR_OUT
cp $PLOTS_PATH $DIR_OUT
cd $DIR_OUT

# Create png files of plots
echo "Extracting plots from archive"
PLOTS_TAR=(${PLOTS_PATH##*/})
if [[ ! -e ${PLOTS_TAR} ]]; then
    echo "ERROR: File ${PLOTS_TAR} does not exist in working directory"
    exit 1
fi
tar -xf $PLOTS_TAR
if [[ `ls -1 *.eps 2>/dev/null | wc -l` -eq 0 ]]; then
    if [[ `ls -1 *.pdf 2>/dev/null | wc -l` -eq 0 ]]; then
	echo "ERROR: No files in .eps or .pdf format found in ${PLOTS_TAR}"
	exit 1
    fi
fi
echo -n "Converting to "
if [[ `ls -1 *.eps 2>/dev/null | wc -l` -ne 0 ]]; then
    echo -n "pdf"
    for eps in `ls -1 *.eps 2>/dev/null`; do
	epstopdf ${eps}
	rm ${eps}
    done
    echo -n " and to "
fi
echo "png format and resizing plots"
for pdf in `ls -1 *${ID_RVs}*.pdf 2>/dev/null`; do
    #pos=(`expr index ${pdf} ${ID_RVs}`)
	new_pdf=`awk 'BEGIN { print substr("'${pdf}'",index("'${pdf}'","'${ID_RVs}'")) }'`
    new_pdf="${DIR_OUT}_${new_pdf}"
    mv ${pdf} ${new_pdf}
    convert ${new_pdf} -resize 400x400 `basename ${new_pdf} .pdf`.png
done
echo "Preparing archive with plots in pdf format"
rm ${PLOTS_TAR}
PLOTS_TAR=(${DIR_OUT}.tar)
tar -cf ${PLOTS_TAR} *.pdf




# ----- Define functions for html output ----------------------

HTML_INDENT=0

function html_newline {
    if [[ $1 -gt 0 ]]; then
	for i in {1..$1} ; do
	    echo "" >> $FILE_OUT
	done
    else
	echo "" >> $FILE_OUT
    fi
}

function html_indentation {
    for i in {1..$HTML_INDENT} ; do
	echo -n " " >> $FILE_OUT
    done
}

function html_line {
    if [[ "$2" == "-" ]]; then
	if [[ $HTML_INDENT -gt 1 ]]; then
	    ((HTML_INDENT=$HTML_INDENT-2))
	fi
    fi
    html_indentation
    echo $1 >> $FILE_OUT

    if [[ "$2" == "+" ]]; then
	((HTML_INDENT=$HTML_INDENT+2))
    fi
}

function html_comment {
    html_line "<!-- $1 -->"
}

function html_tag {
    html_indentation
    echo -n "<$1" >> $FILE_OUT
    if [[ -n ${3} ]]; then echo -n " $3" >> $FILE_OUT; fi
    echo ">$2</$1>" >> $FILE_OUT
}

function html_tag_start {
    html_indentation
    echo -n "<$1" >> $FILE_OUT
    if [[ -n ${3} ]]; then echo -n " $3" >> $FILE_OUT; fi
    echo ">$2" >> $FILE_OUT
    ((HTML_INDENT=$HTML_INDENT+2))
}

function html_tag_end {
    ((HTML_INDENT=$HTML_INDENT-2))
    html_indentation
    echo "</$1>" >> $FILE_OUT
}




# ----- Create webpage ----------------------------------------

# Create html file
echo "Creating html file '${FILE_OUT}'"
touch ${FILE_OUT}
if [[ ! -w ${FILE_OUT} ]]; then
    echo "ERROR: Could not create writable file '${FILE_OUT}' in working directory" 
    exit 1
fi

# Write header and introductory section
html_line "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"
html_newline 2
html_tag_start html
html_tag_start head
html_tag title "${ID} ${JET_ALGO} ${JET_TYPE} jets JEC validation"
html_tag_end head
html_newline 2
html_tag_start body "" "style=\"text-align:left; margin-left:5%; margin-right:5%; font-size:1.2em\""
if [[ ${MODE} == "kalibri" ]]; then
    html_tag h1 "${ID} ${JET_ALGO} ${JET_TYPE} jets Kalibri JEC validation"
else
    html_tag h1 "${ID} ${JET_ALGO} ${JET_TYPE} jets JEC validation"
fi
html_newline
html_comment "Definitions of sample and JEC file etc."
html_line "<p>On this page, validation plots are presented for the ${ID} jet energy corrections for ${JET_ALGO} ${JET_TYPE} jets in CMS. See the top-level <a href=\"https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisJECValidation\">JEC validation</a> TWiki page for a detailed description of the used datasets, the applied event selection, and the validation technique.</p>"
 
# Table listing dataset, JEC tag, etc
html_tag_start table "" "style=\"width:100%; text-align:left; font-size:1em\""

html_tag_start tr
html_tag td "Dataset" "style=\"width:15%\""
html_tag td "${DATASET}" "style=\"font-family:monospace\""
html_tag_end tr

html_tag_start tr
html_tag td "Corrections" "style=\"width:15%\""
html_tag td "${JEC_TAG}" "style=\"font-family:monospace\""
html_tag_end tr

html_tag_start tr
html_tag td "Jets" "style=\"width:15%\""
html_tag td "${JET_ALGO} ${JET_TYPE}" "style=\"font-family:monospace\""
html_tag_end tr

html_tag_start tr
html_tag td "Validation date" "style=\"width:15%\""
html_tag td "`date`" "style=\"font-family:monospace\""
html_tag_end tr

html_tag_start tr
html_tag td "Validation tool" "style=\"width:15%\""
html_tag td "<a href=\"https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration\">Kalibri</a>" "style=\"font-family:monospace\""
html_tag_end tr

html_tag_end table


# Create list of available categories of validation plots

# Check presence of different categories (>0: present, 0: not present)
HAS_PLOTS_ResponseVsGenPt=`ls -1 *${ID_RVsGenPt}*${ID_Response}*.png 2>/dev/null | wc -l`
HAS_PLOTS_ResponseVsEta=`ls -1 *${ID_RVsEta}*${ID_Response}*.png 2>/dev/null | wc -l`
HAS_PLOTS_ResolutionVsGenPt=`ls -1 *${ID_RVsGenPt}*${ID_Resolution}*.png 2>/dev/null | wc -l`
HAS_PLOTS_ResolutionVsEta=`ls -1 *${ID_RVsEta}*${ID_Resolution}*.png 2>/dev/null | wc -l`

# Write table of plots to html file
html_newline 3
html_comment "List of available validation plots"
html_line "<hr>"
html_tag h2 "List of validation plots" "id=\"ListOfValidationPlots\""
html_tag_start ol
if [[ ${HAS_PLOTS_ResponseVsGenPt} -ne 0 ]]; then
    html_tag li "<a href=\"#${ID_RVsGenPt}_${ID_Response}\">${TITLE_ResponseVsGenPt}</a>"
fi
if [[ ${HAS_PLOTS_ResponseVsEta} -ne 0 ]]; then
    html_tag li "<a href=\"#${ID_RVsEta}_${ID_Response}\">${TITLE_ResponseVsEta}</a>"
fi
if [[ ${HAS_PLOTS_ResolutionVsGenPt} -ne 0 ]]; then
    html_tag li "<a href=\"#${ID_RVsGenPt}_${ID_Resolution}\">${TITLE_ResolutionVsGenPt}</a>"
fi
if [[ ${HAS_PLOTS_ResolutionVsEta} -ne 0 ]]; then
    html_tag li "<a href=\"#${ID_RVsEta}_${ID_Resolution}\">${TITLE_ResolutionVsEta}</a>"
fi
html_tag_end ol
html_line "All validation plots are also <a href=\"${PLOTS_TAR}\">archived in pdf format</a>."

# Write plots to html file
function tableOfPlots {
    NUM_UNZOOMED=`ls -1 *${1}*${2}*${4}*[^zoom].png 2>/dev/null | wc -l`;
    NUM_ZOOMED=`ls -1 *${1}*${2}*${4}*[zoom].png 2>/dev/null | wc -l`;
    ((NUM=${NUM_UNZOOMED}+${NUM_ZOOMED}))
    if [[ ${NUM} -ne 0 ]]; then
	COUNT=0
	html_tag_start ol "" "start=\"${5}\""  
	html_tag_start li "<h3 id=\"${1}_${2}\">$3</h3>"
	html_tag_start table "" "style=\"text-align:left\" border=\"0\""
	
	for plot in `ls -1 *${1}*${2}*${4}*.png 2>/dev/null`; do
	    if [[ ${NUM_UNZOOMED} -ne 0 ]]; then
		if [[ ${plot} == *zoom* ]]; then
		    continue
		fi
	    fi
	    html_comment "$4 bin ${COUNT}"
	    ((COUNT=${COUNT}+1))

	    html_tag_start tr
	    
	    html_tag_start td
	    html_line "<img src=\"${plot}\" style=\"border:0\">"
	    html_tag_end td
	    html_tag_start td
	    html_line "( <a href=\"${plot}\">png</a> | <a href=\"`basename ${plot} .png`.pdf\">pdf</a> )"
	    html_tag_end td

	    if [[ ${NUM_ZOOMED} -eq 0 ]]; then
		html_tag_end tr
		continue
	    elif [[ ${NUM_UNZOOMED} -eq 0 ]]; then
		html_tag_end tr
		continue
	    fi
		
	    plot=`basename ${plot} .png`_zoom.png
	    html_tag td "" "style=\"width:5%\""
	    html_tag_start td
	    html_line "<img src=\"${plot}\" style=\"border:0\">"
	    html_tag_end td
	    html_tag_start td
	    html_line "( <a href=\"${plot}\">png</a> | <a href=\"`basename ${plot} .png`.pdf\">pdf</a> )"
	    html_tag_end td
	    	    
	    html_tag_end tr

	done
	html_tag_end table
	html_line "<a href=\"#ListOfValidationPlots\">Back</a> to list of validation plots"
	html_tag_end li
	html_tag_end ol
	html_line "<hr>"
    fi
}

html_newline 3
html_comment "Validation plots in two column format (larger and zoomed y axis range)"
html_line "<hr>"
html_tag h2 "Validation plots"

COUNT_CATEGORIES=0
# Response vs genJetPt
if [[ ${HAS_PLOTS_ResponseVsGenPt} -ne 0 ]]; then
    ((COUNT_CATEGORIES=${COUNT_CATEGORIES}+1))
    html_newline
    tableOfPlots "${ID_RVsGenPt}" "${ID_Response}" "${TITLE_ResponseVsGenPt}" "${ID_Eta}" "${COUNT_CATEGORIES}"
fi
# Response vs eta
if [[ ${HAS_PLOTS_ResponseVsEta} -ne 0 ]]; then
    ((COUNT_CATEGORIES=${COUNT_CATEGORIES}+1))
    html_newline
    tableOfPlots "${ID_RVsEta}" "${ID_Response}" "${TITLE_ResponseVsEta}" "${ID_GenPt}" "${COUNT_CATEGORIES}"
fi
# Resolution vs genJetPt
if [[ ${HAS_PLOTS_ResolutionVsGenPt} -ne 0 ]]; then
    ((COUNT_CATEGORIES=${COUNT_CATEGORIES}+1))
    html_newline
    tableOfPlots "${ID_RVsGenPt}" "${ID_Resolution}" "${TITLE_ResolutionVsGenPt}" "${ID_Eta}" "${COUNT_CATEGORIES}"
fi
# Resolution vs eta
if [[ ${HAS_PLOTS_ResolutionVsEta} -ne 0 ]]; then
    ((COUNT_CATEGORIES=${COUNT_CATEGORIES}+1))
    html_newline
    tableOfPlots "${ID_RVsEta}" "${ID_Resolution}" "${TITLE_ResolutionVsEta}" "${ID_GenPt}" "${COUNT_CATEGORIES}"
fi

html_tag span "<a href=\"https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration\">University of Hamburg calibration group</a>,  `date`" "style=\"font-size:1em\""
html_tag_end body
html_tag_end html

echo "Done"
