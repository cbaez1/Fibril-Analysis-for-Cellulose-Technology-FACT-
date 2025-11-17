(* ::Package:: *)

BeginPackage["FACTpkgFile`"]

impOrig::usage="impOrig[fn_String]\n imports the original images and only selects first channel";
psImportFunc::usage="psImportFunc[path_String]\n Imports the list of pixel sizes for various file types.";
sortFunc1::usage="sortFunc1[fn_]\n Sorts file names based on last number in the name.";
sortFunc2::usage="sortFunc2[fn_]\n Sorts file names based on last number in the name. Uses DigitCharacter..";
dimFunc::usage="dimFunc[data_]\n Outputs the dimensions of either an array or image.";
trimmingFunc::usage="trimmingFunc[images_List, dims_List]\n If the lenghts of the input lists are the same
the function threads ImageTrim.";
pixCountFunc::usage="pixCountFunc[data_List,number_List]\n Counts the number of pixels in 
a scale bar";
tableDisplay1::usage="tableDisplay1[filenames_List]\n prints a table showing the file base names and 
the length of the import file list";
tableDisplay2::usage="tableDisplay2[data_List]\n Prints a table that shows the pixel size associated with 
each image.";
tableDisplay3::usage="tableDisplay3[filenames_, data_List]\n Prints a table showing the number of branch points 
in each skeleton image.";
tableDisplay4::usage="tableDisplay4[filenames_, data_List,ssdata_List]\n Prints a table showing number of junction 
points and total skeleton segment count";
imgDisplay1::usage="imgDisplay1[image_List, size_, ar_]\n dynamic display of images";
formatForExp::usage="formatForExp[data_List, filenames_List]\n formats the width data to export as an excel file";
expToExcel::usage="expToExcel[data_List, filenames_List]\n This version corrects the errors of fromatForExp.";
exportFunc::usage="exportFunc[dir_,name_,data_]\n exports data as one file";
expMult::usage="expMultiple[dir_String, name_List, ext_String, data_List]\n exports data as multiple files";
exportAll::usage="exportAll[datapath_,fnoriginal_,datafn_,kernsize_,mincount_,minelongation_,dropfrac_,
anglethreshold_,ssmincount_,widthdata_,widthmeandat_,refwidthdata_,refwidthmeandat_,expformat_]\n Applies correct export function 
based on expformat input.";
normalizingFunc::usage="normalizingFunc[images_List]\n takes a list of images and normalizes them so that 
all values all lie within 0 and 1.";
thresholdFunc::usage="thresholdFunc[images_List, threshval_List]\n thresholds and binarizes gray images.";
scaleFunc::usage="scaleFunc[data_List,pixsize_,mode_String]\n Multiplies data by the pixel size. Also multiplies
data by 2 if yes was input in mode";
scaleFuncFilt::usage="scaleFuncFilt[data_List,pixsize_,mode_String]\n Filters width data and then multipies by
pixel size."
imgData::usage="imgData[skeleton_Image]\n extracts the skeleton data from the skeleton image.";
imgDataV2::usage="imgDataV2[sscdata_List,binimg_Image]\n Extracts distance transfrom values using the coordinates
from sscdata.";
skeletonFunc::usage="skeletonFunc[binimg_,modbin_]\n Computes the distance transform of the modified binary image and iterates the 
thinning operator on the modified binary until convergence. Then multiplies both results together.";
skeletonFuncV2::usage="skeletonFuncV2[binimg_,modbin_]\n Computes the distance transform of the original binary image and 
iterates the thinning operator on the modified binary until convergence. Then multiplies both results together.";
sklBranchPts::usage="sklBranchPts[skl_List]\n Finds skeleton branch points. Outputs a binary image.";
grayJunctPts::usage="grayJunctPts[skeletonlist_List,branchptsimgs_List,kernsize_]\n Uses the binary image of 
the junction points to find the junction points gray values in the skeleton image.";
branchCompFunc::usage="branchCompFunc[sklimg_]\n Uses ComponentMeasurements to compute the count, perimeter positions,
centroid, medoid, and lable of each branch in a skeleton image.";
posCompile2D::usage="posCompile2D\n Compiled version of Position for 2D images.";
brnchMeasureFunc::usage="brnchMeasureFunc[coords_List, index_List]\n Computes skeleton segment count,centroid, medoid,
and label then formats data in a list.";
nestFunc::usage="nestFunc[data_,x_,checklist_]\n Finds the three nearest points that are within 1 pixel radius of x";
testFunc::usage="testFunc[data_]\n Logical test used to terminate NestWhileList function";
branchPtsSort::usage="branchPtsSort[data_,medoid_]\n Sorts branch points in sequential order.";
dropFunc::usage="dropFunc[data_,dropcount_]\n Used to trim points from the ends of a branch.";
dropFuncV2::usage="dropFuncV2[data_,dropcount_]\n New version of dropFunc where the minimum number
of points per skeleton segment is one instead of zero.";
dropFuncV3::usage="dropFuncV3[data_,dropcount_,ssendpt_,endpttaglen_]\n This version of dropFunc allows
for the user to select whether or not to trim from the side of skeleton segment that is adjacent to 
a skeleton segment end";
ssTrim::usage="ssTrim[data_,dropfrac_,ssendpt_]\n Applies dropFuncV3 and returns a formated result list.";
branchModFunc::usage="branchModFunc[image_,data_,dropfrac_]\n Applies branchPtsSort and dropFunc.";
branchModFuncV2::usage="branchModFuncV2[data_,disttrans_,endpts_,junctpts_]\n Applies branchPtsSort and taggs
each skeleton segment with applicable corresponding junction and end points."
mapModFunc::usage="mapModFunc[data_,dropfrac_]\n Maps branchModFunc.";
mapModFuncV2::usage="mapModFuncV2[ssc_,binimg_,endptsimg_,junctptsimgs_]\n Maps branchModFuncV2";
lineFit3D::usage="lineFit3D[data_]\n formats the data as {x,y,z} where x and y are image coordinates 
and z is the distance transform radial distance then finds the best fit line. Returns the standardized 
branch points, line graphics, vector angle, and branch lable.";
lineFit2D::usage="lineFit2D[data_]\n Formats the data as {n,r} where n is the branch point location 
and r is the distance transform value, then finds the best fit line using singular value decomposition.
Returns the standardized branch points, infinite line graphics, vector angle, and branch lable.";
corrFact::usage="corrFact[angl_]\n Function that computes the skeleton segment orientation corrective factor";
orientCorrect::usage="orientCorrect[data_]\n Determines orientation of skeleton segment then resamples data according 
to the corrective factor";
applyOrientCorr::usage="applyOrientCorr[data_List,usecorr_]\n Applies skeleton segment orientation correction if user selects 
to do so.";
filterJPs::usage="filterJPs[data_,jp1_,jp2_,ep_,joindist_]\n Used to remove spurious junction points. 
These are caused by very small skeleton segments about length 1 to 3 pixels where there are instances when 
subtracting the dilated skeleton junction point (with the intent to separate the skeleton network)
has the uninteded effect of deleting these short skeleton segments. This leaves behind skeleton junction points 
that are appear to not be related to any skeleton segment.";
dragZoomShow::usage="dragZoomShow[Dynamic[rangeVar_],expr_Graphics,opts:OptionsPattern[]]\n Allows the user
to use the mouse to drag and zoom into a graphics image.";
thinTdtV1::usage="thinTdtV1[binimg_,modedbinimg_,binsize_,scalefactor_,maxiter_,stepiter_,size_]\n Outputs the distance
transform image and the corresponding value histogram. Allows you to scroll throught the different thinning iterations
to see the effect on the skeleton image and histogram.";
thinTdtV2::usage="thinTdtV2[img_,binimg_,modedbinimg_,maxiter_,stepiter_]\n Shows the effect of the subsequent 
iterations of the thinning opeartor on the moded binary image. The result is overlayed on top of the active binary image";
thinTdtV3::usage="thinTdtV3[img_,binimg_,modedbinimg_,maxiter_,stepiter_,dimfrac_]\n Used to visualize the effect of the 
thinning operator on partitioned images.";
thinTdtV4::usage="thinTdtV4[img_,binimg_,partmodbins_,maxiter_,dimfrac_,scale_]\n Shows an overlay of the 
skeleton points over the partitioned images. You can hover the mouse over the points to see the width values.";
thinTdtV5::usage="thinTdtV5[img_,binimg_,partmodbins_,maxiter_,dimfrac_,part_,imgsize_:500]\n Use this function 
to export partition images overlaid with the skeleton";
histoDisplay::usage="histoDisplay[data_,fn_,size_,binsize_]\n Displays a histogram of the angle that represents the average rate of change 
of the distance transform values of a branch.";
pltRangeV2::usage="pltRangeV2[data_]\n Changes the plot range values for lineFitDisplay.";
lineFitDisplay::usage="lineFitDisplay[data_,filenames_]\n Outputs a 3D plot showing the branch points and best fit line";
lineFitGraphic::usage="lineFitGraphic[fitres_List,ar_,size_]\n Displays best fit line to branch distance transform data.";
posOfMax::usage="posOfMax[data_,image_]\n Outputs the coordinate of the branch point with the maximum width value";
posOfMean::usage="posOfMean[data_,image_]\n Outputs the coordinate of the branch point closest to the mean branch 
width value";
posOfMedian::usage="posOfMedian[data_,image_]\n Outputs the coordinate of the branch point closest to the median 
branch width value";
nearestCoord::usage="nearestCoord[edgepix_,branchpos_,num_]\n Returns both the input branch position coordinate and the
position of its nearest edge pixel.";
highlightLocations::usage="highlightLocations[branchdat_List,separatebrnch_Image,trimmedimg_Image,binaryimg_Image,loctype_]\n 
Outputs an overlay of the original image, filtered skeleton points, location of point of interets in branch, and a line between
the location of the point of interest and its neareast background pixel.";
highlightLocationsV2::usage="highlightLocationsV2[branchdat_List,separatebrnch_Image,trimmedimg_Image,binaryimg_Image,
loctype_,sklptsize_,locptsize_,circle_String]\n Outputs an overlay of the original image, filtered skeleton points, 
location of point of interets in branch, a line between the location of the point of interest and its neareast background pixel, 
and the user can select to also overlay an enclosing circle with center in point of interest and radius equal to the eucledian disatance
between the point of interest and its nearest background pixel";
scaleBarV1::usage="scaleBar[img_,xfrac_,len_,yfrac_,res_,size_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]\n Adds a 
scale bar to the graphic based on absolute scale bar size in um. (This version does not work)";
scaleBarV2::usage="scaleBarV2[img_,xfrac_,yfrac_,lenfrac_,imgdim_List,pixsize_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]\n 
Adds a scale bar to the graphic based on fraction of scale bar length. (This version has issues)";
scaleBarV3::usage="scaleBarV3[img_,xfrac_,yfrac_,lenfrac_,imgdim_List,pixsize_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]\n 
This version seems to have fixed the scale bar issue.";
resampleFunc::usage="resampleFunc[dat_List,max_Integer]\n This sample is used to reduce the sampling rate of the skeleton segments
before using it in the color maps.";
widthColorMap::usage="widthColorMap[img_,skldat_,scale_,multby2_,colscheme_,max_,units_String,multfactor_,colmaprange_,ptsize_,
legendscale_,lengendnumscale_]\n Displays a raster image of a skeleton segment width color map.";
juncPtGraphic::usage="juncPtGraphic[binimg_Image,junctionpts_List,modbcfilt_,maxrad_,brnchptssize_,junctptsize_,imgsize_,
sklcolor_,JPcolor_]\n Displays a graphic highlighting the skeleton branches and junction points.";
downSamp::usage="downSamp[data_,partlen_]\n This function downsamples the number of points per skeleton segment. Use this 
when the number of points is too large to rasterize. Input an integer greater than 1 to lower the number of points per segment
proportional to that number.";
mapMonitorV1::usage="mapMonitorV1[func_,binimgs_List]\n A version of MapIndexed that outputs a progress bar as it evaluates
the input list.";
mapMonitorV2::usage="mapMonitorV2[func_,skls_List,bc_List,dropfrac_]\n A second version of MapIndexed that outputs a progress 
bar as it evaluates the input list.";
mapMonitorV3::usage="mapMonitorV3[ssc_List,binimgs_List,endptsimg_List,junctptsimgs_List]\n Applies a table that evaluates
mapModFuncV2 and monitors its progress.";
formatDat1::usage="formatDat1[data_,mainlables_,sublables_,subgrouplen_]\n formats the width data to use in 
statistical tables and box-whisker charts.";
outlierFrac::usage="outlierFrac[data_]\n Calculates the fraction of outlier width points.";
statsTable::usage="statsTable[data_Association]\n outputs various statistical values of the data";
statsTablesJoined::usage="statsTablesJoined[data_Association]\n Flatten version of statsTable";
formatDat2::usage="formatDat2[dat1_,dat2_]\n Formats the data to use in statsTableV2";
statsTableV2::usage="statsTableV2[data_Association,grouplengths_]\n This table gathers and presents the statitical data based on 
assign groups";
printStats::usage="printStats[data_,type_,groupnames_,sigdigit_]\n prints the statistics table depending on type selected.";



Begin["`Private`"];


(* ::Section::Closed:: *)
(*Importing images*)


(* ::Subsection::Closed:: *)
(*Import functions*)


(* ::Code::Initialization::Bold:: *)
impOrig[fn_String] := Block[{img},
  img = Import[fn];
  (*selecting only first channel if imported image is multichannel*)
  Which[
  Head[img]=== Image,
  First @ ColorSeparate[img]
  ]
]

psImportFunc[path_String] := Which[
  FileExtension[path] === "txt",
    Import[path, "CSV"]
  ,
  FileExtension[path] === "xlsx",
    Import[path]
  ,
  FileExtension[path] === "mx",
    Import[path]
]

sortFunc1[fn_] := Block[{x},
  ToExpression @ StringCases[FileBaseName[fn], "Cropped_" ~~ x__ -> x
    ]
]

sortFunc2[fn_] := Block[{x},
  ToExpression @ Last @ StringCases[FileBaseName[fn], DigitCharacter..
    ]
]

dimFunc[data_]:=Map[If[Head[#]===Image,ImageDimensions[#],Dimensions[#]]&,data,{1}]



(* ::Subsection::Closed:: *)
(*Trimming function*)


(* ::Code::Initialization::Bold:: *)
trimmingFunc[images_List, dims_List] := 
	If[Length[dims] === Length[images],
    MapThread[ImageTrim[#1, #2]&, {images, dims}],
   {"Lists of unequal length. The lenth of images is " <> ToString[Length[images]] <> 
    " and the length of dimensions is " <> ToString[Length[dims]] <> ". Edit the 
	list of dimensions so that both list are of equal length."}
	]



(* ::Subsection::Closed:: *)
(*Measuring pixel size*)


(* ::Code::Initialization::Bold:: *)
pixCountFunc[data_List,number_List]/;(Length[number]===2):=Block[{fp,sp,fpm,spm},
fp=Select[data,#[[5]]===number[[1]]&][[1,4]];
sp=Select[data,#[[5]]===number[[2]]&][[1,4]];
{fpm,spm}=Sort[{fp,sp},#1[[2,1]]<#2[[2,1]]&];
spm[[2,1]]-fpm[[1,1]]
]

pixCountFunc[data_List,number_List]/;(Length[number]===1||Length[number]===0):=Block[{fp},
fp=Select[data,#[[5]]===number[[1]]&][[1,4]];
fp[[2,1]]-fp[[1,1]]
]



(* ::Subsection::Closed:: *)
(*Display functions 1*)


(* ::Code::Initialization::Bold:: *)
tableDisplay1[filenames_List] := TableForm[Insert[Map[{FileNameTake[#,
   -1]}&, filenames], Length[filenames], {1, 2}], TableHeadings -> {None,
   {"File names", "File count"}}, TableAlignments -> Center]

tableDisplay2[data_List] := TableForm[data, TableHeadings -> {Range[Length[data]], {
  "File names", "Pixel size"}}, TableAlignments -> Center]

tableDisplay3[filenames_, data_List] := TableForm[Thread[{FileBaseName
   /@ filenames, Map[If[Length[#]>0, Keys @ Last[#], 0]&, data]}], TableHeadings ->
   {None, {"File names", "Number of junction points"}}, TableAlignments ->
   Left]
   
   tableDisplay4[filenames_, jp_List, jpref_, ssdata_List] :=  TableForm[Thread[{FileBaseName
   /@ filenames, Map[If[Length[#]>0, Keys @ Last[#], 0]&, jp], Length/@jpref ,Map[Length,ssdata],
   Map[Length@Flatten[#[[All,3]],1]&,ssdata],Map[Total[#[[All,1]]]&,ssdata]}], TableHeadings ->
   {Range[Length[filenames]], {"File names", "# jp original","# jp ref","# Skl segments",
   "Skl counts (Resampled)", "Skl counts"}}, 
   TableAlignments -> Left]

imgDisplay1[image_List, size_, ar_] := If[Head[image[[1]]] === Image,
  Manipulate[Show[image[[i]], ImageSize -> size, AspectRatio-> ar], {i, 1, Length[image
    ], 1}, ControlPlacement -> Top]
  ,image[[1]]
  ]
  


(* ::Subsection::Closed:: *)
(*Export functions*)


(* ::Code::Initialization::Bold:: *)
(*Deprecated. This function does not work as intended.*)
formatForExp[data_List, filenames_List] := Block[{n, paddata, moddat},
   n = Max @ Map[Length, data]; 
   paddata = Map[PadRight[#, n]&, data];
   moddat = Prepend[Transpose[paddata], filenames]; 
   moddat /. {0 -> Nothing}]

(*Use this function for formatting data to export into excel. You can use ctrl + shift + down arrow
 in excel to select all populated cells in a column. Once the data is exported use the first row of index
  positions to sort the columns in the way the order they appear in FACT. For this first select all the 
  data then go to Sort & Filter, then click options and select sort left to right then in Sort by select 
  Row 1 and in order choose smallest to larges, then click Ok.*)
expToExcel[data_List, filenames_List] := Module[{ moddat,moddatindx,datsort,paddata},
moddat=MapThread[Prepend[#1,#2]&,{data,filenames}];
moddatindx=MapThread[Prepend[#1,#2]&,{moddat,Range[Length[filenames]]}];
datsort=Sort[moddatindx,Length[#1]>Length[#2]&];
paddata =PadRight[datsort];
Transpose[paddata]/.{0->Nothing}
]

exportFunc[dir_, name_, data_] := Block[{}, 
	SetDirectory[dir]; 
	Export[name, data]; 
	ResetDirectory[];
	]

expMult[dir_String, name_List, ext_String, data_List] := Block[{}, 
	SetDirectory[dir]; 
	MapThread[Export[FileBaseName[#1] <> ext, #2]&, {name, data}]; 
	ResetDirectory[];
	]
	
exportAll[datapath_,fnoriginal_,datafn_,kernsize_,mincount_,minelongation_,dropfrac_,
anglethreshold_,ssmincount_,widthdata_,widthmeandat_,refwidthdata_,refwidthmeandat_,expformat_]:=Module[{},

Which[
expformat==="xlsx"||expformat==="XLSX",

(*Exporting all skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_allSklWidth."<>
ToString[expformat],expToExcel[widthdata,FileBaseName/@fnoriginal]];
(*Exporting all skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_meanAllSklWidth."<>
ToString[expformat],expToExcel[widthmeandat,FileBaseName/@fnoriginal]];
(*Exporting refined skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_"<>ToString[dropfrac]<>"df_"<>
ToString[anglethreshold]<>"deg"<>"_"<>ToString[ssmincount]<>"smC"<>"_refinedSklWidth."<>
ToString[expformat],expToExcel[refwidthdata,FileBaseName/@fnoriginal]];
(*Exporting mean refined skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_"<>ToString[dropfrac]<>"df_"<>
ToString[anglethreshold]<>"deg"<>"_"<>ToString[ssmincount]<>"smC"<>"_meanRefinedSklWidth."<>
ToString[expformat],expToExcel[refwidthmeandat,FileBaseName/@fnoriginal]],

expformat!="xlsx"||expformat!="XLSX",

(*Exporting all skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_allSklWidth."<>
ToString[expformat],widthdata];
(*Exporting all skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_meanAllSklWidth."<>
ToString[expformat],widthmeandat];
(*Exporting refined skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_"<>ToString[dropfrac]<>"df_"<>
ToString[anglethreshold]<>"deg"<>"_"<>ToString[ssmincount]<>"smC"<>"_refinedSklWidth."<>
ToString[expformat],refwidthdata];
(*Exporting mean refined skeleton width data as xlsx file.*)
exportFunc[First@FileNames["Statistical results",datapath,2],datafn<>"_"<>ToString[kernsize]<>"kS"<>
"_"<>ToString[mincount]<>"mC"<>"_"<>ToString[minelongation]<>"mE"<>"_"<>ToString[dropfrac]<>"df_"<>
ToString[anglethreshold]<>"deg"<>"_"<>ToString[ssmincount]<>"smC"<>"_meanRefinedSklWidth."<>
ToString[expformat],refwidthmeandat]
]
]
  


(* ::Section::Closed:: *)
(*Segmenting images*)


(* ::Subsection::Closed:: *)
(*Normalizing images*)


(* ::Input::Initialization:: *)
normalizingFunc[images_List]:=Block[{imgdata,min,max},
imgdata=Map[If[Head[#]===Image,ImageData[#],#]&,images];
{min,max}=Transpose@Map[MinMax,imgdata];MapThread[Image[#1+#3/(#2-#1)]&,{min,max,imgdata}]
]



(* ::Subsection::Closed:: *)
(*Thresholding*)


(* ::Code::Initialization::Bold:: *)
thresholdFunc[images_List, threshval_List] := 
	If[Length[threshval] === Length[images],
    MapThread[If[#3==="yes",ColorNegate@Binarize[#1, #2],Binarize[#1, #2]]&, 
    {images, threshval[[All,1]], threshval[[All,2]]}],
    "Lists of unequal length. The lenth of normalizedimgs is " <> ToString[Length[images]] <> 
    " and the length of threshvals is " <> ToString[Length[threshval]] <> ". Edit the 
	list of threshvals so that both list are of equal length."
	]
	


(* ::Section::Closed:: *)
(*Image Skeleton*)


(* ::Subsection::Closed:: *)
(*Scaling data*)


(* ::Input::Initialization:: *)
scaleFunc[data_List,pixsize_,mode_String]:=
If[mode==="yes",
Map[If[#>1.,2 pixsize*#,pixsize*#]&,data],data]/.{{}->Nothing}

scaleFuncFilt[data_List,pixsize_,mode_String]:=Module[{select},
select=Select[data,#>=1&];
If[mode==="yes",
Map[If[#>1.,2 pixsize*#,pixsize*#]&,select],select]/.{{}->Nothing}
]

imgData[skeleton_Image]:=Flatten[ImageData[skeleton]]/.{0.:>Nothing,0:>Nothing}

imgDataV2[sscdata_List,binimg_Image]:=Module[{dt,iv},
dt=DistanceTransform[binimg];
iv=Map[ImageValue[dt,#]&,sscdata[[All,2]],{1}];
iv/.{0.:>Nothing,0:>Nothing}
]



(* ::Subsection::Closed:: *)
(*Thinning operation*)


(* ::Input::Initialization:: *)
skeletonFunc[binimg_,modbin_]:=Block[{dt},
dt=DistanceTransform[modbin];
ImageMultiply[dt,Thinning[modbin,Method->"Morphological"]]
]

skeletonFuncV2[binimg_,modbin_]:=Block[{dt},
Thinning[modbin,Method->"Morphological"]
]



(* ::Subsection::Closed:: *)
(*Skeleton segment components*)


(* ::Input::Initialization:: *)
sklBranchPts[skl_List]:=Block[{selcomp,branchptsimg,endptsimg,count},selcomp=Map[SelectComponents[Binarize[#,0.0001],#Count>2&]&,skl];branchptsimg=Map[MorphologicalTransform[#,"SkeletonBranchPoints"]&,selcomp];
endptsimg=Map[MorphologicalTransform[#,"SkeletonEndPoints"]&,selcomp];count=Map[ComponentMeasurements[#,"Count"]&,branchptsimg];{branchptsimg,endptsimg,count}]

grayJunctPts[skeletonlist_List,branchptsimgs_List,kernsize_]:=MapThread[SelectComponents[{#1,MorphologicalComponents[Dilation[#2,kernsize],CornerNeighbors->False]},#Count>1&][[1]]&,{skeletonlist,branchptsimgs}]

(*Trying to replace this slow function. Evaluating PerimeterPositions is very slow. Deprecated*)
(*branchCompFunc[sklimg_]:=Values@ComponentMeasurements[{sklimg,MorphologicalComponents[Binarize[sklimg,0.0001]]},{"Count","PerimeterPositions","Centroid","Medoid","Label"},#Count>0&,CornerNeighbors->True]*)

posCompile2D=Compile[{{array, _Integer, 2}, {indx, _Integer, 0}}, 
 Position[array, indx], Parallelization -> True, CompilationTarget -> "C", 
 RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}];

(*This function replaced branchCompFunc. Its much faster.*)
  brnchMeasureFunc[coords_List, index_List] /;(Length[coords[[1,1]]]===2):= 
  Module[{count, centroid,medoid, lable, datlist}, 
  count = Length /@ coords; 
  centroid = Map[{Total[#[[All, 1]]], Total[#[[All, 2]]]} / Length[#]&, coords]; 
 medoid = First /@ MapThread[Nearest[#1, #2, DistanceFunction -> Automatic]&, {coords, centroid}]; 
  datlist = Thread[{count, coords, centroid, medoid, index}]; 
  {#1, {#2}, #3, #4, #5}& @@@ datlist
  ]



(* ::Subsection::Closed:: *)
(*Skeleton segment trimming and filtering*)


(* ::Input::Initialization:: *)
nestFunc[data_,x_,checklist_]:=Block[{xmod,checklistmod,near,filt,checklnew},
xmod=If[Depth[x]>2,x[[1]],x];
checklistmod=If[Depth[x]>2,x[[2]],checklist];
near=Nearest[data,xmod,{3,1},DistanceFunction->ChessboardDistance];
filt=Flatten@DeleteElements[near,checklistmod];
checklnew=AppendTo[checklistmod,filt];
{filt,checklnew}
]

testFunc[data_]:=UnsameQ[{},data[[1]]]

branchPtsSort[data_,medoid_]:=Block[{checklist,near,path1start,path1res,path2start,path2res,path1mod,path2mod},
checklist={medoid};
near=Nearest[data,medoid,3];
path1start=near[[2]];
AppendTo[checklist,path1start];
path1res=NestWhileList[nestFunc[data,#1,checklist]&,path1start,testFunc[#]&];
path2start=near[[3]];
AppendTo[checklist,path2start];
path2res=NestWhileList[nestFunc[data,#1,checklist]&,path2start,testFunc[#]&];
path1mod=Join[{path1res[[1]]},Rest[path1res[[All,1]]]]/.{}->Nothing;
path2mod=Join[{path2res[[1]]},Rest[path2res[[All,1]]]]/.{}->Nothing;
Which[
Length[data]===3,
{near[[2]],near[[1]],near[[3]]},
Length[data]===4,
Union@Join[Reverse[path2mod],{medoid},path1mod],
(*This corrects for when the medoid of the branch happens to be an end point.*)
Length[data]>4&&Rest[path1mod]!=path2mod,
Join[Reverse[path2mod],{medoid},path1mod],
Length[data]>4&&Rest[path1mod]===path2mod,
Join[{medoid},path1mod]
]
]

dropFunc[data_,dropcount_]:=Drop[Drop[data,dropcount[[1]]],-dropcount[[2]]]

(*In this version the minimum amount of points per skeleton segment is 1. The previous version allowed the minimum to be zero which was not my intention.*)
dropFuncV2[data_,dropcount_]:=Module[{d1},
d1=Drop[Drop[data,dropcount[[1]]],-dropcount[[2]]];
If[Length[d1]===0,Drop[Drop[data,dropcount[[1]]],-dropcount[[2]]+1],d1]
]

(*Currently this function only removes half of the drop fraction for the skeleton segment side associated with an end (when ssendpt = 1). Not sure if this is what is intended. Alternatively, when ssendpt= 1 I could double the drop count for skeleton segments containing an end point.*)
dropFuncV3[data_,dropcount_,ssendpt_,endpttaglen_]:=Module[{d1},
d1=If[ssendpt===2,
Drop[Drop[data,dropcount[[1]]],-dropcount[[2]]],
Which[
endpttaglen[[1]]===0&&endpttaglen[[2]]===0,
Drop[Drop[data,dropcount[[1]]],-dropcount[[2]]],
endpttaglen[[1]]===2&&endpttaglen[[2]]===0,
Drop[data,-dropcount[[2]]],
endpttaglen[[1]]===0&&endpttaglen[[2]]===2,
Drop[data,dropcount[[1]]],
endpttaglen[[1]]===2&&endpttaglen[[2]]===2,
data
]
];
If[Length[d1]===0,
{Median[data]},
d1]
]

ssTrim[data_,dropfrac_,ssendpt_]:=Module[{dropcount,trimskl,trimdtvals},
dropcount={Floor[#],Ceiling[#]}&[data[[1]]*dropfrac/2];
trimskl=dropFuncV3[data[[2]],dropcount,ssendpt,Length/@data[[7]]];
trimdtvals=dropFuncV3[data[[3]],dropcount,ssendpt,Length/@data[[7]]];
{Length[trimskl],trimskl,trimdtvals,data[[4]],data[[5]],data[[6]],data[[7]],data[[8]]}
]

(*Deprecated*)
branchModFunc[image_,data_,dropfrac_]:=Block[{pixpos,distsort,imgval,res1,dropcount,trimskl,trimdtval},
pixpos=Union[Flatten[data[[2]],1]];
distsort=If[Length[pixpos]>2,branchPtsSort[pixpos,data[[4]]],pixpos];
imgval=ImageValue[image,distsort];
res1={data[[1]],distsort,imgval,data[[3]],data[[4]],data[[5]]};
dropcount={Floor[#],Ceiling[#]}&[res1[[1]]*dropfrac/2];
trimskl=dropFuncV2[distsort,dropcount];
trimdtval=dropFuncV2[imgval,dropcount];
{Length[trimskl],trimskl,trimdtval,data[[3]],data[[4]],data[[5]]}
]

branchModFuncV2[data_,disttrans_,endpts_,junctpts_]:=Module[{sspos,sssort,distvals,endpttag,junctpttag},
(*Skeleton segment positions.*)
sspos=Flatten[data[[2]],1];
(*Skeleton segment sorted*)
sssort=If[Length[sspos]>2,branchPtsSort[sspos,data[[4]]],sspos];
(*Distance transform values*)
distvals=ImageValue[disttrans,sssort];
(*Tagging which skeleton segment end has an end-point associated with it*)
endpttag=Flatten/@Map[Nearest[endpts,#,{1,2}]&,{First[#],Last[#]}&@sssort];
(*Tagging which skeleton segment end has an junction point associated with it. The max distance an untrimmed skeleton segment end can be from a junction point is Sqrt[2^2+2^2] = 2.83*)
junctpttag=If[
junctpts==={},
{{},{}},
Nearest[junctpts,{First[#],Last[#]}&@sssort,{All,2.85}]
];
(*Formatting data*)
{Length[sspos],sssort,distvals,data[[3]],data[[4]],data[[5]],endpttag,junctpttag}
]

(*Deprecated*)
mapModFunc[data_,dropfrac_]:=Map[branchModFunc[data[[1]],#,dropfrac]&,data[[2]]]

mapModFuncV2[ssc_,binimg_,endptsimg_,junctptsimgs_]:=Module[{distimg,endpts,junctpts},
{distimg,endpts,junctpts}={DistanceTransform[binimg],ImageValuePositions[endptsimg,1],ImageValuePositions[junctptsimgs,1]};
Map[branchModFuncV2[#,distimg,endpts,junctpts]&,ssc]
]

lineFit3D[data_]:=Block[{normalizedpts,linepts,line,vector,vectangl},
normalizedpts=Standardize[#,Mean,1&]&@Map[Flatten,Thread[{data[[2]],data[[3]]}]];
linepts={Mean[normalizedpts],Flatten[Last[SingularValueDecomposition[normalizedpts,1]]]};
line=InfiniteLine@@linepts;
vector=linepts[[1]]-linepts[[2]];
vectangl=VectorAngle[vector,ReplacePart[vector,3->0]]*180/Pi;
{normalizedpts,line,vectangl,data[[6]]}
]

lineFit2D[data_]:=Module[{normalizedpts,linepts,line,vector,vectangl},
normalizedpts=Standardize[#,Mean,1&]&@Map[Flatten,Thread[{Range[Length[data[[3]]]],data[[3]]}]];
linepts={Mean[normalizedpts],Flatten[Last[SingularValueDecomposition[normalizedpts,1]]]};
line=InfiniteLine@@linepts;
vector=linepts[[1]]-linepts[[2]];
vectangl=VectorAngle[vector,ReplacePart[vector,2->0]]*180/Pi;
{normalizedpts,line,vectangl,data[[6]]}
]

corrFact[angl_]:=If[angl<45,Cos[angl Degree],Sin[angl Degree]];

orientCorrect[data_]:=Block[{standardized,coords,line,vector,vecangl,gph2d,anglandcountRR0norm,corrfact,resampval,resampdat},
standardized=Standardize[data[[2]],Mean,1&];
coords={Mean[standardized],Flatten[Last[SingularValueDecomposition[standardized,1]]]};
line=InfiniteLine@@coords;
vector=coords[[1]]-coords[[2]];
vecangl=(VectorAngle[vector,ReplacePart[vector,2->0]]*180/Pi)/.Indeterminate->0;
corrfact=corrFact[vecangl];
resampval=Ceiling[data[[1]]/corrfact];
resampdat=ArrayResample[data[[3]],resampval];
ReplacePart[data,3->resampdat]
]

applyOrientCorr[data_List,usecorr_]/;(usecorr==="yes"):=Map[If[#[[1]]>1,orientCorrect[#],#]&,data]

applyOrientCorr[data_List,usecorr_]/;(usecorr==="no"):=data



(* ::Subsection::Closed:: *)
(*Junction point refining*)


(* ::Input::Initialization:: *)
filterJPs[data_,jp1_,jp2_,ep_,joindist_]:=Module[{search1,remove1,newjp1,search2,remove2,newjp2},
(*Searching for junction points that have at least two skeleton segment ends within a radius of 2.85. Use this to remove the junction points with only 1 skeleton end segment withing a radius of 2.85. Intended to remove spurious junction points.*)
search1=Union[ReplaceAll[Nearest[Flatten[Map[{First[#],Last[#]}&,data[[All,2]]],1],jp1,{All,2.85}],{}->Nothing]];
(*Selecting from search the junction points that only have 1 skeleton segment end withing a radius of 2.85 pix*)
remove1=Flatten[Select[search1,Length[#]===1&],1];
(*Remove the junction point from jskls that are nearest to the remove1 skeleton segment end list.*)
newjp1=Select[jp2,FreeQ[Flatten[Nearest[jp2,remove1,{1,2.85}],1],#]&];
(*Searches for the junction points that have two end points within a 2.85 pixel radius. Intended to remove another type of spurious junction points.*)
search2=Union[ReplaceAll[Nearest[ep,jp1,{All,2.85}],{}->Nothing]];
remove2=Map[First,Select[search2,Length[#]>1&]];
newjp2=Select[jp2,FreeQ[Flatten[Nearest[jp2,remove2,{1,2.85}],1],#]&];
Sort[Map[Mean,Gather[newjp2,Norm[#1-#2]<joindist&]],#1[[1]]<#2[[1]]&]
]



(* ::Subsection::Closed:: *)
(*Display functions 2*)


(* ::Input::Initialization:: *)
dragZoomShow[Dynamic[rangeVar_],expr_Graphics,opts:OptionsPattern[]]:=
DynamicModule[{ start, end,image}, DynamicWrapper[EventHandler[Dynamic[Show[start; image, Epilog -> If[ValueQ[start], {Opacity[0], EdgeForm[Directive[Orange, Dashed,Thickness[.005]]], Rectangle[start, MousePosition["Graphics"]]}, {}],opts]], {"MouseDown" :> (start = MousePosition["Graphics"]), "MouseUp" :> (end = MousePosition["Graphics"]; If[MatchQ[start, {_Real, _Real}] && MatchQ[end, {_Real, _Real}], 
rangeVar={Sort[{start[[1]], end[[1]]}],Sort[{start[[2]], end[[2]]}]}];
 If[start[[1]] === end[[1]] || start[[2]] === end[[2]], rangeVar=PlotRange[expr]]; start =. )}],

image = If[MatchQ[rangeVar,{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}],
Show[expr, PlotRange->rangeVar,opts],(rangeVar=PlotRange[#];#)&@expr
],
TrackedSymbols:>{rangeVar}]];

dragZoomShow[expr_,opts:OptionsPattern[]]:=
DynamicModule[{rangeVar},dragZoomShow[Dynamic[rangeVar],expr,opts]];

thinTdtV1[binimg_,modedbinimg_,binsize_,scalefactor_,maxiter_,stepiter_,size_]:=DynamicModule[{dt,min,max,fold,mult,hist},
dt=DistanceTransform[binimg];
{min,max}=MinMax[dt];
fold=FoldList[Thinning[#1,#2]&,modedbinimg,Table[stepiter,{maxiter/stepiter}]];
mult=ImageMultiply[#,dt]&/@fold;
hist=Table[Histogram[(Map[If[#>1.,2scalefactor #,scalefactor #]&,Flatten[ImageData[mult[[i]]]]])/.{0.:>Nothing},{binsize},"Probability",ImageSize->IntegerPart[size*0.6],PlotRange->All,AxesLabel->{"Width","Probability"}],{i,1,Length[mult],1}];
Manipulate[
GraphicsRow[{
Show[ImageAdjust[min+mult[[i]]/(max-min),{0,0,0.5}],ImageSize->size*1],
hist[[i]]},10,ImageSize->size],
{i,1,Length[mult],1},ControlPlacement->Top]
]

thinTdtV2[img_,binimg_,modedbinimg_,maxiter_,stepiter_]:=DynamicModule[{fold},
fold=FoldList[Thinning[#1,#2]&,modedbinimg,Table[stepiter,{maxiter/stepiter}]];
Manipulate[
Show[HighlightImage[img,{EdgeForm[],Red,fold[[i]]}],ImageSize->size],
{i,1,Length[fold],1},{{size,500},300,1000,100},ControlPlacement->Top]
]

thinTdtV2Video[img_,binimg_,modedbinimg_,maxiter_,stepiter_]:=DynamicModule[{dt,min,max,fold,mult},
dt=DistanceTransform[binimg];
{min,max}=MinMax[dt];
fold=FoldList[Thinning[#1,#2]&,modedbinimg,Table[stepiter,{maxiter/stepiter}]];
mult=Binarize[ImageMultiply[#,min+dt/(max-min)],0.001]&/@fold;
Information[Video[#,Appearance->"Minimal"]&@Animate[
Show[HighlightImage[img,{EdgeForm[{Thick,Red}],mult[[i]]}],ImageSize->500],
{i,1,Length[mult],1},ControlType->None,Paneled->False],"VideoTracks"]
]

thinTdtV2VideoBin[img_,binimg_,modedbinimg_,maxiter_,stepiter_]:=DynamicModule[{dt,min,max,fold,mult},
dt=DistanceTransform[binimg];
{min,max}=MinMax[dt];
fold=FoldList[Thinning[#1,#2]&,modedbinimg,Table[stepiter,{maxiter/stepiter}]];
mult=Binarize[ImageMultiply[#,min+dt/(max-min)],0.001]&/@fold;
Information[Video[#,Appearance->"Minimal"]&@Animate[
Show[HighlightImage[binimg,{EdgeForm[{Thick,Red}],mult[[i]]}],ImageSize->500],
{i,1,Length[mult],1},ControlType->None,Paneled->False],"VideoTracks"]
]

thinTdtV2VideoJoined[img_,binimg_,modedbinimg_,maxiter_,stepiter_]:=DynamicModule[{dt,min,max,fold,mult},
dt=DistanceTransform[binimg];
{min,max}=MinMax[dt];
fold=FoldList[Thinning[#1,#2]&,modedbinimg,Table[stepiter,{maxiter/stepiter}]];
mult=Binarize[ImageMultiply[#,min+dt/(max-min)],0.001]&/@fold;
Information[Video[#,Appearance->"Minimal"]&@Animate[
GraphicsRow[{
Show[HighlightImage[img,{EdgeForm[{Thick,Red}],mult[[i]]}],ImageSize->500],
Show[HighlightImage[binimg,{EdgeForm[{Thick,Red}],mult[[i]]}],ImageSize->500]}],
{i,1,Length[mult],1},ControlType->None,Paneled->False],"VideoTracks"]
]

thinTdtV3[img_,binimg_,modedbinimg_,maxiter_,stepiter_,dimfrac_]:=DynamicModule[{dims,partimgs,partbins,partmodbins,dt,min,max,mult},
dims=ImageDimensions[img];
partimgs=Flatten@ImagePartition[img,dims/dimfrac];
partbins=Flatten@ImagePartition[binimg,dims/dimfrac];
partmodbins=Flatten@ImagePartition[modedbinimg,dims/dimfrac];
dt=DistanceTransform/@partbins;
{min,max}=MinMax[Flatten[ImageData/@dt]];
Manipulate[
mult=ImageMultiply[Thinning[partmodbins[[i]],k],min+dt[[i]]/(max-min)];
Show[HighlightImage[If[showbinimg==True,partbins[[i]],partimgs[[i]]],{EdgeForm[{Thin,Red}],Opacity[0.2],If[delsmall==True,DeleteSmallComponents[#],#]&@Pruning[Binarize[mult,0.01],p]}],ImageSize->size],
{i,1,Length[partbins],1},{k,1,maxiter,stepiter},{p,0,20,5},{{size,500},300,1000,100},{showbinimg,{False,True}},{delsmall,{False,True}},ControlPlacement->Top]
]

thinTdtV4[img_,binimg_,partmodbins_,maxiter_,dimfrac_,scale_]:=DynamicModule[{dims,partimgs,partbins,partmodbinimgs,dt,skeleton,pixpos,diameter,radial,partsklt},
dims=ImageDimensions[img];
dt=DistanceTransform[binimg];
skeleton=ImageMultiply[Thinning[partmodbins,maxiter],dt];
partsklt=Flatten@ImagePartition[skeleton,dims/dimfrac];
partimgs=Flatten@ImagePartition[img,dims/dimfrac];
partbins=Flatten@ImagePartition[binimg,dims/dimfrac];
partmodbinimgs=Flatten@ImagePartition[partmodbins,dims/dimfrac];
pixpos=Map[PixelValuePositions[Binarize[#,0.00001],1]&,partsklt];
diameter=Map[NumberForm[#,{3,2}]&,Table[Map[If[#>1.,2scale #,scale #]&, PixelValue[partsklt[[i]],pixpos[[i]]]],{i,Length[partsklt]}],{2}];
radial=Map[NumberForm[#,{3,2}]&,Table[scale PixelValue[partsklt[[i]],pixpos[[i]]],{i,Length[partsklt]}],{2}];
Manipulate[
HighlightImage[
Which[
imgtype==1,
ImageAdjust@partimgs[[i]],
imgtype==2,
partbins[[i]],
imgtype==3,
partmodbinimgs[[i]]
],{Red,PointSize[Small],MapThread[Tooltip[#1,#2]&,{Point/@(pixpos[[i]]-0.5),If[rad==True,radial[[i]],diameter[[i]]]}]},ImageSize->size],
{i,1,Length[partbins],1},{{size,500},300,1000,100},{imgtype,{1,2,3}},{rad,{False,True}},ControlPlacement->Top]
]

thinTdtV5[img_,binimg_,partmodbins_,maxiter_,dimfrac_,part_,imgsize_:500]:=
Block[{dims,partimgs,partbins,dt,min,max,skeleton,pixpos,partsklt},
dims=ImageDimensions[img];
dt=DistanceTransform[binimg];
{min,max}=MinMax[dt];
skeleton=ImageMultiply[Thinning[partmodbins,maxiter],dt];
partsklt=Flatten@ImagePartition[skeleton,dims/dimfrac];
partimgs=Flatten@ImagePartition[img,dims/dimfrac];
partbins=Flatten@ImagePartition[binimg,dims/dimfrac];
pixpos=Map[PixelValuePositions[Binarize[#,0.00001],1]&,partsklt];
{HighlightImage[partbins[[part]],{Red,PointSize[Small],Point[pixpos[[part]]-0.5]},ImageSize->imgsize],
Show[ImageAdjust@partimgs[[part]],Graphics[{Red,PointSize[Small],Point[pixpos[[part]]]}],ImageSize->imgsize]}
]

histoDisplay[data_,fn_,size_,binsize_]:=
Histogram[data[[All,3]],{binsize},"Probability",AxesLabel->{"Angle (deg)","Probability"},PlotRange->{{0,55},All},PlotLabel->FileBaseName[fn],ImageSize->size]

pltRangeV2[data_]:=Block[{xmin,xmax,ymin,ymax,zmin,zmax,delta},
{xmin,xmax}=MinMax@data[[All,1]];
{ymin,ymax}=MinMax@data[[All,2]];
{zmin,zmax}=MinMax@data[[All,3]];
delta={xmax-xmin,ymax-ymin,zmax-zmin};
{{xmin,xmin+Max[delta]},{ymin,ymin+Max[delta]},{zmin,zmin+Max[delta]}}
]

lineFitDisplay[data_,filenames_]:=
Graphics3D[{{Directive[AbsolutePointSize[6],Black],Point[data[[1]]]},{Directive[AbsoluteThickness[4],Red],data[[2]]}},Axes->True,PlotRange->pltRangeV2[data[[1]]](*Automatic*),AxesLabel->{"x","y","z"},AspectRatio->Automatic,PlotLabel->ToString[FileBaseName[filenames]]<>"\nBranch lable: "<>ToString[data[[4]]]<>" -> "<>ToString[Round[data[[3]],0.01]]<>" deg"
]

lineFitGraphic[fitres_List,ar_,size_]:= Graphics[{Directive[AbsolutePointSize[8], Black],
 Point[fitres[[1]]], Directive[AbsoluteThickness[3], Red], fitres[[2]]},
 Axes -> True, PlotRange -> All, AspectRatio -> ar, 
 PlotLabel -> ToString[fitres[[4]]] <> " --> " <> ToString[Round[fitres[[3]], 0.1]],ImageSize->size
 ]

posOfMax[data_,image_]:=Block[{coords,pv,mix},
coords=data[[2]];
pv=PixelValue[image,coords];
mix=Sort[Thread[{pv,coords}],#1[[1]]>#2[[1]]&];
First[mix][[2]]
]

posOfMean[data_,image_]:=Block[{coords,pv,mean,nearestindex},
coords=data[[2]];
pv=PixelValue[image,coords];
mean=Mean[pv];
nearestindex=First@Nearest[pv->"Index",mean,1];
coords[[nearestindex]]
]

posOfMedian[data_,image_]:=Block[{coords,pv,median,nearestindex},
coords=data[[2]];
pv=PixelValue[image,coords];
median=Median[pv];
nearestindex=First@Nearest[pv->"Index",median,1];
coords[[nearestindex]]
]

nearestCoord[edgepix_,branchpos_,num_]:=Block[{nearest},
nearest=Map[Nearest[edgepix,#,num,DistanceFunction->EuclideanDistance]&,branchpos];
Flatten[Table[Map[{branchpos[[i]],#1}&,nearest[[i]]],{i,Length[nearest]}],1]
]

highlightLocations[branchdat_List,separatebrnch_Image,trimmedimg_Image,binaryimg_Image,loctype_,sklptsize_,locptsize_]:=Block[{maxpos,meanpos,medianpos,edgepts,ncmax,ncmean,ncmedian,ncmedoid},
maxpos=Map[posOfMax[#,separatebrnch]&,branchdat];
meanpos=Map[posOfMean[#,separatebrnch]&,branchdat];
medianpos=Map[posOfMedian[#,separatebrnch]&,branchdat];
edgepts=ImageValuePositions[EdgeDetect[binaryimg],1];
ncmax=nearestCoord[edgepts,maxpos,1];
ncmean=nearestCoord[edgepts,meanpos,1];
ncmedian=nearestCoord[edgepts,medianpos,1];
ncmedoid=nearestCoord[edgepts,branchdat[[All,5]],1];
HighlightImage[trimmedimg,{Green,PointSize[sklptsize],Point[Flatten[branchdat[[All,2]],1]],
Red,Thickness[0.001],
Which[loctype==="max",
Line/@ncmax,
loctype==="mean",
Line/@ncmean,
loctype==="median",
Line/@ncmedian,
loctype==="medoid",
Line/@ncmedoid],
Blue,PointSize[locptsize],
Which[loctype==="max",
Point/@maxpos,
loctype==="mean",
Point/@meanpos,
loctype==="median",
Point/@medianpos,
loctype==="medoid",
Point/@branchdat[[All,5]]]},ImageSize->500]
]

highlightLocationsV2[branchdat_List,separatebrnch_Image,trimmedimg_Image,binaryimg_Image,loctype_,sklptsize_,locptsize_,circle_String]:=
Block[{pos,edgepts,nc,rad},
pos=Which[
loctype==="max",
Map[posOfMax[#,separatebrnch]&,branchdat],
loctype==="mean",
Map[posOfMean[#,separatebrnch]&,branchdat],
loctype==="median",
Map[posOfMedian[#,separatebrnch]&,branchdat],
loctype==="medoid",
branchdat[[All,5]]];
edgepts=ImageValuePositions[EdgeDetect[binaryimg],1];
nc=nearestCoord[edgepts,pos,1];
rad=EuclideanDistance[#1,#2]&@@@nc;
HighlightImage[trimmedimg,{Green,PointSize[sklptsize],Point[Flatten[branchdat[[All,2]],1]],Thickness[0.001],
Blue,PointSize[locptsize],
Point/@pos,
If[circle==="yes",Thread[Circle[pos,rad]]],Red,Line/@nc},ImageSize->500]
]

scaleBarV1[img_,xfrac_,yfrac_,len_,imgdim_List,pixsize_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]:=Show[img,Graphics[{color,Thickness[linethick],Line[{{xfrac*imgdim[[1]]-(len/pixsize)/2,yfrac*imgdim[[2]]},{xfrac*imgdim[[1]],yfrac*imgdim[[2]]},{xfrac*imgdim[[1]]+(len/pixsize)/2,yfrac*imgdim[[2]]}}],Style[Text[StringJoin[ToString[len],unit],{xfrac*imgdim[[1]],yfrac*imgdim[[2]]-ynumpos}],fontsize,Bold]}]]

(*Deprecated*)
scaleBarV2[img_,xfrac_,yfrac_,lenfrac_,imgdim_List,pixsize_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]:=Block[{newlen,xmin,xmax,ypos,bkgoffset,graphic},
newlen=IntegerPart[lenfrac*imgdim[[1]]*pixsize];
{xmin,xmax}={xfrac*imgdim[[1]]-(newlen/pixsize)/2,xfrac*imgdim[[1]]+(newlen/pixsize)/2};
ypos=yfrac*imgdim[[2]];
bkgoffset=IntegerPart[lenfrac*imgdim[[1]]*0.02];
(*Not sure why these offsets are necessary. Here I am computing the offset as a fraction of the image column size. Here we first plot a white background for the scalebar, then the scalebar, then the size text.*)
graphic=Graphics[{White,Thickness[linethick+0.01],
Line[{{xmin-bkgoffset+0.003*imgdim[[1]],ypos},{xmax+bkgoffset-0.01*imgdim[[1]],ypos}}],color,Thickness[linethick],Line[{{xmin+0.003*imgdim[[1]],ypos},{xmax-0.01*imgdim[[1]],ypos}}],
Style[Text[StringJoin[ToString[NumberForm[newlen,{6,2}]],unit],{xfrac*imgdim[[1]],ypos-ynumpos}],fontsize,Bold,Background->White]},ImageSize->imgdim];ImageResize[HighlightImage[img,graphic],imgdim]
]

scaleBarV3[img_,xfrac_,yfrac_,lenfrac_,imgdim_List,pixsize_,color_,ynumpos_,fontsize_Integer,linethick_,unit_String]:=Module[{len,xmin,xmax,ypos,bkgoffset,graphic1,newimg,newxmin,newxmax,newlen,printlength,graphic2},
len=Round[lenfrac*imgdim[[1]]*pixsize,.1];
{xmin,xmax}={xfrac*imgdim[[1]]-(len/pixsize)/2,xfrac*imgdim[[1]]+(len/pixsize)/2};
ypos=yfrac*imgdim[[2]];
bkgoffset=IntegerPart[lenfrac*imgdim[[1]]*0.02];
graphic1=Graphics[{White,Thickness[linethick],Line[{{xmin,ypos},{xmax,ypos}}]}];
newimg=Binarize@ImageResize[HighlightImage[ConstantImage[0,imgdim],graphic1],imgdim];
{newxmin,newxmax}=First[Values[ComponentMeasurements[newimg,"BoundingBox"]]][[All,1]];
newlen=(newxmax-newxmin)*pixsize;
printlength=If[newlen<1,Round[newlen,0.01],Round[newlen,1]];
graphic2=Graphics[{White,Thickness[linethick+0.01],
Line[{{xmin-bkgoffset,ypos},{xmax+bkgoffset,ypos}}],color,Thickness[linethick],
Line[{{xmin,ypos},{xmax,ypos}}],
Style[Text[StringJoin[ToString[printlength],unit],{xmin+(xmax-xmin)/2,ypos-ynumpos}],fontsize,Bold,Background->White]}];
ImageResize[HighlightImage[img,graphic2],imgdim]
(*{len,newlen,len/newlen//N}*)
]

resampleFunc[dat_List,max_Integer]:=Module[{part},
part=Partition[dat,UpTo[max]];
Map[Extract[#,Ceiling[Length[#]/2]]&,part]
]

widthColorMap[img_,skldat_,scale_,multby2_,colscheme_,max_,units_String,multfactor_,colmaprange_,ptsize_,legendscale_,lengendnumscale_]:=Module[{pixpos,vals,valsscaled,range},
pixpos=Flatten[#,1]&@Map[resampleFunc[#,max]&,skldat[[All,2]]];
vals=Flatten[#,1]&@Map[resampleFunc[#,max]&,skldat[[All,3]]];
valsscaled=multfactor*scaleFunc[vals,scale,multby2];
range=If[Length[colmaprange]>0,colmaprange,{#1,#2}&@@MinMax[valsscaled]];
Rasterize[
Legended[
HighlightImage[img,{PointSize[ptsize],MapThread[{#1,#2}&,{ColorData[{colscheme,range}]/@valsscaled,Point/@pixpos}]},ImageSize->Full],
BarLegend[{colscheme,range},LegendLabel->"Width ("<>units<>")",
LegendMarkerSize->(*Last[ImageDimensions[img]]**)legendscale,LabelStyle->Directive[Black,FontSize->Round[(*Last[ImageDimensions[img]]**)lengendnumscale]]]
],RasterSize->1500]
]

juncPtGraphic[binimg_Image,junctionpts_List,modbcfilt_,maxrad_,brnchptssize_,junctptsize_,imgsize_,sklcolor_,JPcolor_]:=
HighlightImage[binimg,{PointSize[brnchptssize],sklcolor,modbcfilt[[All,2]],JPcolor,PointSize[junctptsize],Point[junctionpts]},ImageSize->imgsize]

downSamp[data_,partlen_]:=Module[{f,l,moddat,partdat},
{f,l}={First[#],Last[#]}&@data;
moddat=If[Length[data]>2,Take[data,{2,-2}],data];
partdat=First/@Partition[moddat,partlen];
Union@Join[{f},partdat,{l}]
]



(* ::Subsection::Closed:: *)
(*Map monitor*)


(* ::Input::Initialization:: *)
mapMonitorV1[func_,binimgs_List]:=Block[{x=0,l=Length[binimgs]},
Monitor[
MapIndexed[(x=First[#2];func[#1[[1]],#1[[2]]])&,binimgs],
{ProgressIndicator[x/l],StringJoin[ToString[x], StringJoin["/",ToString[l]]]}]
]

mapMonitorV2[func_,skls_List,bc_List,dropfrac_]:=Block[{x=0,l=Length[skls]},
Monitor[
MapIndexed[(x=First[#2];func[#1,dropfrac])&,Thread[{skls,bc}]],
{ProgressIndicator[x/l],StringJoin[ToString[x], StringJoin["/",ToString[l]]]}]
]

mapMonitorV3[ssc_List,binimgs_List,endptsimg_List,junctptsimgs_List]:=Module[{x=0,l=Length[binimgs]},
Monitor[
Table[mapModFuncV2[ssc[[x]],binimgs[[x]],endptsimg[[x]],junctptsimgs[[x]]],{x,Length[binimgs]}],
{ProgressIndicator[x/l],StringJoin[ToString[x], StringJoin["/",ToString[l]]]}]
]



(* ::Section::Closed:: *)
(*Width data statistics*)


(* ::Subsection::Closed:: *)
(*Statistical values*)


(* ::Input::Initialization:: *)
formatDat1[data_,mainlables_,sublables_,subgrouplen_]:=Association[
#1->Association[#2]&@@@
Thread[{mainlables,TakeList[MapThread[ToString[#2]->#1&,{data,sublables}],subgrouplen]}]
]

outlierFrac[data_]:=Block[{len1,quart,inbetweendata,len2},
len1=Length/@data;
quart=Quartiles/@data;
inbetweendata=Table[Select[data[[i]],Between[#,{quart[[i,1]]-1.5(quart[[i,3]]-quart[[i,1]]),quart[[i,3]]+1.5(quart[[i,3]]-quart[[i,1]])}]&],{i,Length[data]}];
len2=Length/@inbetweendata;
(1-len2/len1)*100.
]

statsTable[data_Association]:=Block[{statmetrics,outfrac},
statmetrics=Map[Flatten[{Mean[#],Skewness[#],Quartiles[#],StandardDeviation[#],Length[#]}]&,data,{2}];
outfrac=outlierFrac[Flatten[Values[Values[data]],1]];
Flatten/@Transpose[{Flatten[Values[Values[statmetrics]],1],outfrac}]
]

statsTablesJoined[data_Association]:=Block[{statmetrics,outfrac},
statmetrics=Map[Flatten[{Mean[#],Skewness[#],Quartiles[#],StandardDeviation[#],Length[#]}]&,data,{1}];
outfrac=outlierFrac[Values[data]];
Flatten/@Transpose[{Values[statmetrics],outfrac}]
]

formatDat2[dat1_,dat2_]:=MapThread[Append[#1,#2]&,{dat1,dat2}]

statsTableV2[data_Association,grouplengths_]:=Block[{statmetrics,outfrac,list1,modlist1},statmetrics=Map[Flatten[{Mean[#],Skewness[#],Quartiles[#],StandardDeviation[#],Length[#]}]&,data,{2}];outfrac=TakeList[outlierFrac[Flatten[Values[Values[data]],1]],grouplengths];
list1=Apply[Flatten@List[##]&,Values[Map[Normal,statmetrics]],{2}];
modlist1=MapThread[formatDat2[#1,#2]&,{list1,outfrac},1];
{Keys[statmetrics],modlist1}
]

printStats[data_,type_,groupnames_,sigdigit_]:=If[type==="no",
TableForm[Map[NumberForm[#,{4,sigdigit}]&,Transpose/@data[[2]],{3}],TableDirections->Column,TableAlignments->{Center,Center},TableHeadings->{data[[1]],{"Images","Mean","Skewness","1st Quartile","Median","3rd Quartile","StdDev","Count","Outliers(%)"}}],TableForm[Map[NumberForm[#,{3,2}]&,data,{2}],
TableAlignments->Center,TableHeadings->{groupnames,{"Mean","Skewness","1st Quartile","Median","3rd Quartile","StdDev","Count","Outliers(%)"}}]
]



(* ::Section::Closed:: *)
(*End*)


End[];
EndPackage[]
