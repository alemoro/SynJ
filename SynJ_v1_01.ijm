// SynJ - fluoresce synapse analysis in ImageJ

/*

Start developing 2018.10.21
Start release 2019.11.30

Modify
	
*/

var majVer = 1;
var minVer = 01;
var about = "Developed by Alessandro Moro<br>"
			+ "<i>Department of Functional Genomics</i> (FGA)<br>"
			+ "<i>Centre of neuroscience and cognitive research</i> (CNCR)<br>"
			+ "<i>Vrij Universiteit</i> (VU) Amsterdam.<br>"
			+ "<i>email: a.moro@vu.nl</i><br><br><br>";

// Initialize the variables
var mC = 2; // morphology channel
var sC = 3; // synapse channel
var xC = 1; // protein X channel
var minSoma = 6; // minimum soma radius
var maxSoma = 10; // maximim soma radius
var oneSoma = false;
var pxSat = 0.1; // percentage of saturated pixel for the Enhance Contrast
var gbS = 1; // sigma value for Gaussian Blur filter
var lineW = 5; // line width for Ridge detection
var minLength = 20; // minimum neurite length in um (I think)
var highC = 230; // high contrast for Ridge detection
var lowC = 87; // low contrast for Ridge detection
var wth = 6; // radius for the white-tophat filter
var autoT = 6; // radius for the auto local threshold
var synS = 0.5; // sigmal for the Gaussian Blur for synapses
var nThr = 1; // number of stander deviations above the average intensity for synapses detection
var neuritePad = 1; // increase width in neurite width to detect synapses
var nClosing = 1; // number of closing iteration after Ridge detection
var bSomaDetect = false; // detect synapses in the soma?
var bFolder = true; // save RoiManager and data in the same folder as the image
var saveName = "fullname"; // which name to use to save
var morRoiColor = "white"; // color for the outline of the morphology
var keepSynColor = "magenta"; // color for keeping the synapses
var delSynColor = "red"; // color for the synapses to delete
var saveDir = "saveDir";
var nameToSave = "nameToSave";


// Toolset buttons

// leave one empty slot
macro "Unused Tool -1-" {} 

// Soma Analysis -> Automatic detection of soma? Define soma mask
var sCmds1 = newMenu("Detect soma Menu Tool", newArray("Auto detection", "Define soma", "Edit soma", "-", "Get from SynD", "Batch analysis", "-", "Convert *.nd2"));
macro "Detect soma Menu Tool - Ca00 T0b10S T6b10o Tcb10m"{
	readOptions();
	if (nImages == 1) {
		Stack.setDisplayMode("color");
		Stack.setChannel(mC);
	}
	cmd1 = getArgument();
	if (cmd1 == "Auto detection"){
		//waitForUser("Sorry, no automatic detection yet.\nBut this will start the overall analysis");
		autoDetectSoma();
		editSoma();
		roiManager("deselect");
		roiManager("show all without labels");
		run("Select None");
		autoDetectNeurite();
		editNeurites();
		roiManager("select", 1);
		Roi.setStrokeColor(morRoiColor);
		roiManager("update");
		roiManager("deselect");
		roiManager("show all without labels");
		run("Select None");
		autoDetectSynapses();
		deleteWrongSynapses();
		if (roiManager("count") > 2) {
			for (r = 0; r < roiManager("count")-2; r++) {
				roiManager("select", r+2);
				roiManager("rename", "Synapse_"+r+1);
			}
		}
		roiManager("show all without labels");
	}
	else if (cmd1 == "Define soma"){
		setTool("polygon");
		waitForUser("Draw soma");
		roiManager("add");
		roiManager("select", 0);
		roiManager("rename", "Soma");
		roiManager("select", 0);
		Roi.setStrokeColor(morRoiColor);
		roiManager("update");
		roiManager("deselect");
		roiManager("show all without labels");
		run("Select None");
		editSoma();
		roiManager("deselect");
		roiManager("show all without labels");
		run("Select None");
	}
	else if (cmd1 == "Edit Soma"){
		// draw soma for all cells
		editSoma();
	}
	else if (cmd1 == "Get from SynD") {
		imgDir = getDirectory("Select images folder");
		mskDir = getDirectory("Select neuritemask folder");
		imgFile = getFileList(imgDir);
		mskFile = getFileList(mskDir);
		nImg = imgFile.length;
		nMsk = mskFile.length;
		for (i = 0; i < nImg; i++) { // loop throught the images
			if (endsWith(imgFile[i], ".tif")) {
				//open(imgDir + imgFile[i]);
				//title = getTitle();
				//imgID = getImageID();
				for (m = 0; m < nMsk; m++) { // loop to find the right mask
					bImg = startsWith(mskFile[m], imgFile[i]);
					bMsk = indexOf(mskFile[m], "-neuritemask") > 0;
					if (bImg && bMsk) {
						open(imgDir + imgFile[i]);
						title = getTitle();
						imgID = getImageID();
						setBatchMode(true);
						open(mskDir + mskFile[m]);
						mskTitle = getTitle();
						run("Split Channels");
						selectWindow(mskTitle + " (green)"); close(); // remove the confocal component
						selectWindow(mskTitle + " (blue)"); // get the soma
						setOption("BlackBackground", true);
						run("Convert to Mask");
						run("Dilate"); run("Dilate"); // twice for better cleaning
						run("Create Selection");
						roiManager("Add");
						roiManager("Select", 0);
						roiManager("Rename", "Soma");
						Roi.setStrokeColor(morRoiColor);
						selectWindow(mskTitle + " (red)");
						imageCalculator("XOR", mskTitle + " (red)",mskTitle + " (blue)"); // get the neuriteMask
						run("Convert to Mask");
						run("Create Selection");
						roiManager("Add");
						roiManager("Select", 1);
						roiManager("Rename", "Neurites");
						Roi.setStrokeColor(morRoiColor);
						selectWindow(mskTitle + " (red)"); close();
						selectWindow(mskTitle + " (blue)"); close();
						setBatchMode(false);
						autoDetectSynapses(); // detect synapses
						getSaveName(); // and save all
						roiManager("Save", saveDir + "RoiSet_" + nameToSave + ".zip");
						saveIntensity(saveDir, nameToSave);
						roiManager("reset"); // clean the workspace
						selectImage(imgID); close(); // close the image if there is not a mask or when it's done
						if (isOpen("Results")) {
							selectWindow("Results"); run("Close");
						}
					}
				}
				
			}
		}
	}
	else if (cmd1 == "Batch analysis") {
		// first ask for the folder
		imgDir = getDirectory("Select images folder");
		imgFile = getFileList(imgDir);
		nImg = imgFile.length;
		for (i = 0; i < nImg; i++) { // loop throught the images
			if (endsWith(imgFile[i], ".tif")) {
				open(imgDir + imgFile[i]);
				title = getTitle();
				imgID = getImageID();
				// draw soma for all cells
				setTool("polygon");
				waitForUser("Draw soma");
				roiManager("add");
				roiManager("select", 0);
				roiManager("rename", "Soma");
				roiManager("select", 0);
				Roi.setStrokeColor(morRoiColor);
				roiManager("update");
				roiManager("deselect");
				roiManager("show all without labels");
				run("Select None");
				// save the soma before continue
				getSaveName();
				roiManager("Save", saveDir + "RoiSet_" + nameToSave + ".zip");
				// close this image and continue
				selectImage(imgID); close();
				roiManager("reset");
			}
		}
		// now enter in batch mode
		setBatchMode(true);
		for (i = 0; i < nImg; i++) { // loop throught the images
			if (endsWith(imgFile[i], ".tif")) {
				open(imgDir + imgFile[i]);
				title = getTitle();
				imgID = getImageID();
				// open the RoiManager and detect neurites
				roiManager("open", imgDir + "RoiSet_" + replace(title, ".tif", ".zip"));
				autoDetectNeurite();
				roiManager("select", 1);
				Roi.setStrokeColor(morRoiColor);
				roiManager("update");
				roiManager("deselect");
				roiManager("show all without labels");
				run("Select None");
				// and synapses
				autoDetectSynapses();
				deleteWrongSynapses();
				if (roiManager("count") > 2) {
					for (r = 0; r < roiManager("count")-2; r++) {
						roiManager("select", r+2);
						roiManager("rename", "Synapse_"+r+1);
					}
				}
				roiManager("show all without labels");
				// then save everything
				getSaveName();
				roiManager("Save", saveDir + "RoiSet_" + nameToSave + ".zip");
				saveIntensity(saveDir, nameToSave);
				// close image and reset RoiManager
				selectImage(imgID); close();
				roiManager("reset");
			}
		}
		waitForUser("Batch analysis completed");
	} else if (cmd1 == "Convert *.nd2") {
		// Get the main Directory assuming there is a folder per sequence
	workDir = getDirectory("Select main diretory");
	workFold = getFileList(workDir);
	nFold = workFold.length;

	// ask for a saving string
	Dialog.create("Save nd2 sequence");
	Dialog.addString("Condition: ", "date_condition");
	Dialog.addCheckbox("Preserve colors", true);
	Dialog.addCheckbox("Save stack", true);
	Dialog.addCheckbox("Save MAX projection", true);
	Dialog.show;
	saveTitle = Dialog.getString();
	bColor = Dialog.getCheckbox();
	bStack = Dialog.getCheckbox();
	bMAX = Dialog.getCheckbox();
	
	// considering the name of the folder with information for coverslip
	setBatchMode(true);
	for(d=0; d<nFold; d++){
		showStatus("Converting nd2 sequences " + d + "/" + nFold);
		showProgress(d/nFold);
		//waitBar(d, nFold);
		workF = workDir + workFold[d];
		workFiles = getFileList(workF);
		nFile = workFiles .length;
		for(f=0; f<nFile; f++){
			tempFile = workF + workFiles[f];
			if(endsWith(tempFile, "nd2")){
			
				// open nd2 sequences with Bio-Format, all images at one
				run("Bio-Formats Importer", "open=" + tempFile + " autoscale color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

				// get information for saving (save in the same folder)
				savePath = getDirectory("image");
			
				// considering the nd2 sequences are usually "Seq####.nd2 - Seq####.nd2 (series n)
				seq = IJ.pad(f,4); // to have the sequence number
				n = 1; // series number
				nImg = nImages; // total number of images
				while(n <= nImg){
					nn = n;
					if(nImg >= 10) nn = IJ.pad(nn,2); // adjust the series number
					// selectWindow(workFold[d] + "Seq"+seq+".nd2 - Seq"+seq+".nd2 (series "+nn+")");
					selectImage(1);
					rename(saveTitle + "_cs" + (d + 1) + "_" + IJ.pad(n,3));
					tempTitle = getTitle(); // for calling after aligment

					// register stack
					getDimensions(width, height, channels, slices, frames);
					if(channels > 1){
						getDimensions(ww,hh,cc,sc,fc);
						run("Split Channels");
						for(c=1;c<=cc;c++){
							selectWindow("C"+c+"-"+tempTitle);
							//run("StackReg", "transformation=Affine");
						}
						if(cc==3){
							if(bColor){
								run("Merge Channels...", "c1=C2-"+tempTitle+" c2=C1-"+tempTitle+" c3=C3-"+tempTitle+" create");
								rename(tempTitle);
							}else{
								run("Merge Channels...", "c1=C1-"+tempTitle+" c2=C2-"+tempTitle+" c3=C3-"+tempTitle+" create");
							}
						}else if (cc == 4) {
							run("Merge Channels...", "c1=C1-"+tempTitle+" c2=C2-"+tempTitle+" c3=C3-"+tempTitle+" c4=C4-"+tempTitle+" create");
						} else {
							run("Merge Channels...", "c1=C1-"+tempTitle+" c3=C2-"+tempTitle+" create");
						}
						run("Z Project...", "projection=[Max Intensity]");
						rename("MAX_" + tempTitle);
						selectWindow(tempTitle);
						nTemp = 2;
					} else {
						if(slices > 1){
							run("StackReg", "transformation=Affine");
							run("Z Project...", "projection=[Max Intensity]");
							rename("MAX_" + tempTitle);
							selectWindow(tempTitle);
							nTemp = 2;
						} else {
							nTemp = 1;
						}
					}
	
					// save file(s) checking if there is a MAX projection
					if(bStack){
						run("Save", "save=" + savePath + tempTitle + ".tif");
					}
					close();
					if(nTemp == 2){
						selectWindow("MAX_" + tempTitle);
						if(bMAX){
							run("Save", "save=" + savePath + "MAX_" + tempTitle + ".tif");
						}
						close();
					}
					n ++;
				}
				while(nImages>0) close();
			}
		}
	}
	showStatus("Converting nd2 sequences: DONE!");
	return;
	}
	if (nImages == 1) {
		roiManager("select", 0);
		Roi.setStrokeColor(morRoiColor);
		roiManager("update");
		roiManager("deselect");
		roiManager("show all without labels");
		run("Select None");
	}
}

// Neurite Analysis -> Automatic detection of neurites, Define neurites mask
var sCmds2 = newMenu("Detect neurites Menu Tool", newArray("Auto detection", "Edit neurites"));
macro "Detect neurites Menu Tool - C0a0 T0b10N T7b10e Tcb10u"{
	readOptions();
	Stack.setDisplayMode("color");
	Stack.setChannel(mC);
	cmd2 = getArgument();
	if (cmd2 == "Auto detection"){
		autoDetectNeurite();
		editNeurites();
	} else {
		editNeurites();
	}
	roiManager("select", 1);
	Roi.setStrokeColor(morRoiColor);
	roiManager("update");
	roiManager("deselect");
	roiManager("show all without labels");
	run("Select None");
}

// Neurite Analysis -> Automatic detection of neurites, Define neurites mask
var sCmds3 = newMenu("Detect synapses Menu Tool", newArray("Auto detection", "Delete synapses", "-", "Delete all synapses"));
macro "Detect synapses Menu Tool - C00a T0b10S T7b10y Tcb10n"{
	readOptions();
	Stack.setDisplayMode("color");
	Stack.setChannel(sC);
	getPixelSize(unit, pixelWidth, pixelHeight);
	nIter = round(neuritePad/pixelWidth);
	cmd3 = getArgument();
	if (cmd3 == "Auto detection"){
		autoDetectSynapses();
	} else if (cmd3 == "Delete synapses") {
		deleteWrongSynapses();
	} else {
		deleteAllSynapses();
	}
	if (roiManager("count") > 2) {
		for (r = 0; r < roiManager("count")-2; r++) {
			roiManager("select", r+2);
			roiManager("rename", "Synapse_"+r+1);
		}
	}
	roiManager("show all without labels");
}

macro "Mark syanpses Tool - Cf0f R0055R7255R3a55 Cf11 L33ccL43bcLc33cLb33b" {
	setSlice(sC);
	getCursorLoc(x, y, z, flags);
	if (flags == 16){
		nROI = roiManager("count");
		r = 2;
		bROI = 0;
		while((r < nROI) && (bROI == 0)){
			roiManager("Select", r);
			bROI = selectionContains(x, y);
			r++;
		}
		if (bROI == 1) {
			if (matches(Roi.getStrokeColor(), delSynColor)){
				roiManager("Set Color", keepSynColor);
			} else {
				roiManager("Set Color", delSynColor);
			}
		}
	}
}

// Save
var sCmds4 = newMenu("Save Menu Tool", newArray("Save all", "Roi Manager", "Intensity", "Adjust table", "-", "Sholl Analysis", "-", "Collapse data"));
macro "Save Menu Tool - C555D11D12D13D14D15D16D17D18D19D1aD1bD1cD21D27D28D29D2aD2bD2cD2dD31D33D35D37D38D39D3aD3bD3cD3dD3eD41D43D45D47D48D4eD51D53D55D57D58D5aD5bD5cD5eD61D63D65D67D68D6aD6bD6cD6eD71D73D75D77D78D7eD81D83D85D87D88D8eD91D93D95D97D98D9eDa1Da3Da5Da7Da8DaeDb1Db7Db8DbeDc1Dc2Dc3Dc4Dc5Dc6Dc7Dc8Dc9DcaDcbDccDcdDce"{
	readOptions();
	cmd4 = getArgument();
	if (cmd4 == "Roi Manager") {
		getSaveName();
		roiManager("Save", saveDir + "RoiSet_" + nameToSave + ".zip");
	} else if (cmd4 == "Intensity") {
		getSaveName();
		saveIntensity(saveDir, nameToSave);
	} else if (cmd4 == "Save all") {
		getSaveName();
		roiManager("Save", saveDir + "RoiSet_" + nameToSave + ".zip");
		saveIntensity(saveDir, nameToSave);
	} else if (cmd4 == "Adjust table") {
		bCont = true;
		while (bCont) {
			adjustTable();
			bCont = getBoolean("Would you like to adjust another set?");
		}
		return;
	} else if (cmd4 == "Sholl Analysis") {
		// first ask a couple of questions
		Dialog.create("Sholl intensity options");
		Dialog.addChoice("Region to analyze", newArray("Neurite", "Synapses", "All"));
		Dialog.addNumber("Sholl radiud (in um)", 5);
		Dialog.addCheckbox("Batch Analysis", false);
		Dialog.show();
		regSholl = Dialog.getChoice();
		sRadius = Dialog.getNumber();
		bBatch = Dialog.getCheckbox();
		if (bBatch) {
			imgDir = getDirectory("Select images folder");
			imgFile = getFileList(imgDir);
			nImg = imgFile.length;
			setBatchMode(true);
			// count how many images
			nTif = 0;
			for (i = 0; i < nImg; i++) {
				if (endsWith(imgFile[i], ".zip")) {
					nTif++;
				}
			}
			cImg = 0;
			toc = 1;
			for (i = 0; i < nImg; i++) { // loop throught the images
				if (endsWith(imgFile[i], ".tif")) {
					cImg++;
					if (File.exists(imgDir + "RoiSet_" + replace(imgFile[i],".tif",".zip"))) {
						tic = getTime();
						print("Image "+cImg+"/"+nTif+" estimate "+toc*(nTif-cImg+1)+" minutes");
						open(imgDir + imgFile[i]);
						title = getTitle();
						imgID = getImageID();
						// Open the Roi Manager
						roiManager("open", imgDir + "RoiSet_" + replace(imgFile[i],".tif",".zip"));
						getSaveName();
						SynJ_Sholl(regSholl, sRadius);
						saveAs("Results_"+nameToSave, saveDir + "/Sholl_" + nameToSave + ".csv");
						selectWindow("Results"); run("Close");
						roiManager("reset");
						selectImage(imgID); close();
						toc = (getTime() - tic)/1000/60;
					}
				}
			}
		} else {
			SynJ_Sholl(regSholl, sRadius);
		}
	} else if (cmd4 == "Collapse data") {
		workDir = getDirectory("Choose a Directory");
		workFiles = getFileList(workDir);
		nFile = workFiles.length;
		// get the number of result files
		nCSV = 0;
		for (i = 0; i < nFile; i++) {
			if (endsWith(workFiles[i],"csv")) {
				nCSV++;
			}
		}
		// assign the new variables
		allTitle = newArray(nCSV); //Cell ID
		allCsID = newArray(nCSV); // coverslip ID
		allExpID = newArray(nCSV); // experiment ID, condition, protein
		allGrpID = newArray(nCSV);
		allSomaA = newArray(nCSV); // Soma area
		allSomaM = newArray(nCSV); // somatic morphology intensity
		allSomaS = newArray(nCSV); // somatic synaptic marker
		allSomaX = newArray(nCSV); // somatic other marker
		allNeuriteL = newArray(nCSV); // Neurite length
		allNeuriteE = newArray(nCSV); // Neurite ends
		allNeuriteM = newArray(nCSV); // neurite morphology
		allNeuriteS = newArray(nCSV); // neurite synaptic marker
		allNeuriteX = newArray(nCSV); // neurite other marker
		allSynapseA = newArray(nCSV); // mean synapses area
		allSynapseM = newArray(nCSV); // mean intensity synaptic morphology
		allSynapseS = newArray(nCSV); // mean intensity synaptic marker
		allSynapseX = newArray(nCSV); // mean intensity synaptic other marker
		allSynapseD = newArray(nCSV); // synapses per um
		// collect the data
		csv=0;
		for (i = 0; i < nFile; i++) {
			if (endsWith(workFiles[i],"csv")) {
				open(workDir+workFiles[i]);
				//Table.rename(workFiles[i], "Results");
				meanM = 0;
				meanS = 0;
				meanX = 0;
				meanA = 0;
				nRes = getValue("results.count");
				for (j = 2; j < nRes; j++) {
					meanA = meanA + getResult("Area", j);
					meanM = meanM + getResult("Morphology", j);
					meanS = meanS + getResult("Synapses", j);
					meanX = meanX + getResult("Other", j);
				}
				//print(nRes);
				// considering data saved as *[something_]*csID_cellID
				nameParts = split(workFiles[i], "_");
				csN = nameParts.length-2; // -1 to start the array at 0, additional -1 to get the second to last part
				//nameParts[3] = substring(nameParts[3], 0, lengthOf(nameParts[3])-4);
				nameParts[5] = substring(nameParts[5], 0, lengthOf(nameParts[5])-4);
				tempTitle = nameParts[2]+"_"+nameParts[3]+"_"+nameParts[4]+"_"+nameParts[5];
				// get the number of ends from the brackets in the neurite region
				nEnd = getResultString("Region", 1);
				if (lengthOf(nEnd)>10) {
					nEnd = substring(nEnd, 10, lengthOf(nEnd)-1);
				} else {
					nEnd = 0;
				}
				// store the data
				//allTitle[csv] = substring(workFiles[i], 0, lengthOf(workFiles[i])-4);
				allTitle[csv] = tempTitle;
				allCsID[csv] = nameParts[csN];
				allExpID[csv] = nameParts[3];
				allGrpID[csv] = substring(nameParts[csN], 2);
				allSomaA[csv] = getResult("Area", 0);
				allSomaM[csv] = getResult("Morphology", 0);
				allSomaS[csv] = getResult("Synapses", 0);
				allSomaX[csv] = getResult("Other", 0);
				allNeuriteL[csv] = getResult("Area", 1);
				allNeuriteE[csv] = nEnd;
				allNeuriteM[csv] = getResult("Morphology", 1);
				allNeuriteS[csv] = getResult("Synapses", 1);
				allNeuriteX[csv] = getResult("Other", 1);
				allSynapseA[csv] = meanA / (nRes-2);
				allSynapseM[csv] = meanM / (nRes-2);
				allSynapseS[csv] = meanS / (nRes-2);
				allSynapseX[csv] = meanX / (nRes-2);
				allSynapseD[csv] = (nRes-2) / getResult("Area", 1);
				csv++;
				run("Close");
			}
		}
		
		//selectWindow("Results");run("Close");
		// save new data
		for (i = 0; i < nCSV; i++) {
			setResult("Cell ID", i, allTitle[i]);
			setResult("Coverslip ID", i, allCsID[i]);
			setResult("Experiment ID", i, allExpID[i]);
			setResult("Group ID", i, parseInt(allGrpID[i]));
			setResult("Soma area", i, allSomaA[i]);
			setResult("Soma morphology", i, allSomaM[i]);
			setResult("Soma synaptic", i, allSomaS[i]);
			setResult("Soma other", i, allSomaX[i]);
			setResult("Neurite length", i, allNeuriteL[i]);
			setResult("Neurite ends", i, allNeuriteE[i]);
			setResult("Neurite morphology", i, allNeuriteM[i]);
			setResult("Neurite synaptic", i, allNeuriteS[i]);
			setResult("Neurite other", i, allNeuriteX[i]);
			setResult("Synaptic area", i, allSynapseA[i]);
			setResult("Synaptic morphology", i, allSynapseM[i]);
			setResult("Synaptic synaptic", i, allSynapseS[i]);
			setResult("Synaptic other", i, allSynapseX[i]);
			setResult("Synaptic per um", i, allSynapseD[i]);
		}
		updateResults();
	}
}

/*
// Quick plot function
macro "Plot Action Tool - Cf0f" {
	regLabel = newArray("Soma", "Neurite", "Synapses");
	parLabel = newArray("Area", "Morphology", "Synapses", "Other");
	readOptions();
	Dialog.create("Plot SynJ data");
	Dialog.addChoice("Region to plot", regLabel);
	Dialog.addChoice("Parameter to plor", parLabel);
	Dialog.addCheckbox("Normalize?", false);
	Dialog.show();
	regPlot = Dialog.getChoice();
	parPlot = Dialog.getChoice();
	bNorm = Dialog.getCheckbox();
	// Use the loaded collapse data
	
}
*/

// Options
macro "Options... Action Tool - C77bD3eD4eD5eD6bD6cD6dD7aD89D98Da7Db6Dc6Dd6De4De5D2aD5dDa2Dd5D59D68D69D77D78D86D87D96D1aD1bD1cD29D2bD39D49D4bD4cD4dD58D67D76D85D92D93D94Da1Db1Db2Db4Dc1Dc4Dd4De3D5aD6aD79D88D95D97Da5Da6D19D91D4aD5bDa4Db5D3aD5cDa3Dc5"{
	setOptions();
}

// Documentation!!!
/*
macro "Help... Action Tool - C313Dce"{
	message = "<html>"
	 + "<h2>Version " + majVer + "." + minVer + "</h2><br>"
	 + about + "<br>"
	 + "The documentation could be found "
	 + "<a href=\"https://github.com/alemoro/DCVpHluorin_toolseet\">here.</a><br>"
	 + "<a href=\"http://www.johanneshjorth.se/SynD/SynD.html\">SynD</a><br>"
	 + "<a href=\"https://dreamdealer.net/redmine/projects\">Fusion Analysis 2</a><br>"
	 + "<a href=\"http://bigwww.epfl.ch/thevenaz/stackreg/\">StackReg</a><br>"
	 + "<a href=\"http://imagej.net/MorphoLibJ\">MorphoLibJ</a><br>"
	 + "<a href=\"https://imagej.net/Spots_colocalization_(ComDet)\">ComDet</a>.<br>"
	Dialog.create("Help");
	Dialog.addMessage("Version " + majVer + "." + minVer + ", \nclick \"Help\" for more");
	Dialog.addHelp(message);
	Dialog.show;
	//showMessage("Not very usefull", message);
}
*/

// General function
function readOptions() {
	mC = call("ij.Prefs.get", "SynJ.MorphologyChannel", true);
	sC = call("ij.Prefs.get", "SynJ.SynapseChannel", true);
	xC = call("ij.Prefs.get", "SynJ.OtherChannel", true);
	minSoma = call("ij.Prefs.get", "SynJ.minSoma", true);
	maxSoma = call("ij.Prefs.get", "SynJ.maxSoma", true);
	oneSoma = call("ij.Prefs.get", "SynJ.oneSoma", true);
	pxSat = call("ij.Prefs.get", "SynJ.SaturatedPixels", true);
	gbS = call("ij.Prefs.get", "SynJ.GaussinBlurSigma", true);
	lineW = call("ij.Prefs.get", "SynJ.RidgeDetectionLineWidth", true);
	minLength = call("ij.Prefs.get", "SynJ.MinimumLength", true);
	highC = call("ij.Prefs.get", "SynJ.RidgeDetectionHighContrast", true);
	lowC = call("ij.Prefs.get", "SynJ.RidgeDetectionLowContrast", true);
	wth = call("ij.Prefs.get", "SynJ.WhiteTopHatRadius", true);
	autoT = call("ij.Prefs.get", "SynJ.LocalThresholdRadius", true);
	synS = call("ij.Prefs.get", "SynJ.synapseGaussianBlur", true);
	nThr = call("ij.Prefs.get", "SynJ.synapseFilterThreshold", true);
	neuritePad = call("ij.Prefs.get", "SynJ.neuritePadding", true);
	nClosing = call("ij.Prefs.get", "SynJ.closingIteration", true);
	bSomaDetect = call("ij.Prefs.get", "SynJ.detectInSoma", true);
	bFolder = call("ij.Prefs.get", "SynJ.saveInImageFolder", true);
	saveName = call("ij.Prefs.get", "SynJ.saveName", true);
	morRoiColor = call("ij.Prefs.get", "SynJ.morphologyColor", true);
	keepSynColor = call("ij.Prefs.get", "SynJ.keepSynapsesColor", true);
	delSynColor = call("ij.Prefs.get", "SynJ.deleteSynapsesColor", true);
}

function writeOptions() {
	call("ij.Prefs.set", "SynJ.MorphologyChannel", mC);
	call("ij.Prefs.set", "SynJ.SynapseChannel", sC);
	call("ij.Prefs.set", "SynJ.OtherChannel", xC);
	call("ij.Prefs.set", "SynJ.minSoma", minSoma);
	call("ij.Prefs.set", "SynJ.maxSoma", maxSoma);
	call("ij.Prefs.set", "SynJ.oneSoma", oneSoma);
	call("ij.Prefs.set", "SynJ.SaturatedPixels", pxSat);
	call("ij.Prefs.set", "SynJ.GaussinBlurSigma", gbS);
	call("ij.Prefs.set", "SynJ.RidgeDetectionLineWidth", lineW);
	call("ij.Prefs.set", "SynJ.MinimumLength", minLength);
	call("ij.Prefs.set", "SynJ.RidgeDetectionHighContrast", highC);
	call("ij.Prefs.set", "SynJ.RidgeDetectionLowContrast", lowC);
	call("ij.Prefs.set", "SynJ.WhiteTopHatRadius", wth);
	call("ij.Prefs.set", "SynJ.LocalThresholdRadius", autoT);
	call("ij.Prefs.set", "SynJ.synapseGaussianBlur", synS);
	call("ij.Prefs.set", "SynJ.synapseFilterThreshold", nThr);
	call("ij.Prefs.set", "SynJ.neuritePadding", neuritePad);
	call("ij.Prefs.set", "SynJ.closingIteration", nClosing);
	call("ij.Prefs.set", "SynJ.detectInSoma", bSomaDetect);
	call("ij.Prefs.set", "SynJ.saveInImageFolder", bFolder);
	call("ij.Prefs.set", "SynJ.saveName", saveName);
	call("ij.Prefs.set", "SynJ.morphologyColor", morRoiColor);
	call("ij.Prefs.set", "SynJ.keepSynapsesColor", keepSynColor);
	call("ij.Prefs.set", "SynJ.deleteSynapsesColor", delSynColor);
}

function setOptions() {
	// if first time don't read the option
	bFirst = call("ij.Prefs.get", "SynJ.firstRun", true);
	if (bFirst) {
		call("ij.Prefs.set", "SynJ.firstRun", false);
	} else {
		readOptions();
	}
	Dialog.create("SynJ Options");
	Dialog.addMessage("Channels options");
	Dialog.addNumber("Morphology channel", mC);
	Dialog.addNumber("Synapse channel", sC);
	Dialog.addNumber("Other channel", xC);
	Dialog.addMessage("Soma Detection");
	Dialog.addNumber("Minimum soma radius", minSoma);
	Dialog.addNumber("Maximum soma radius", maxSoma);
	Dialog.addCheckbox("Detect one soma", oneSoma);
	Dialog.addMessage("Ridge Detection");
	Dialog.addNumber("% of saturated pixel", pxSat);
	Dialog.addNumber("Gaussian Blur sigma", gbS);
	Dialog.addNumber("Estimate line width", lineW);
	Dialog.addNumber("Minimum neurite length", minLength);
	Dialog.addNumber("High contrast", highC);
	Dialog.addNumber("Low contrast", lowC);
	Dialog.addMessage("Synapse detection");
	Dialog.addNumber("Tophat radius", wth);
	Dialog.addNumber("Synapse radius", autoT);
	Dialog.addNumber("Synapse smooth (0.1-1)", synS);
	Dialog.addNumber("Synapse threshold (std)", nThr);
	Dialog.addNumber("Neurite padding (um)", neuritePad);
	Dialog.addNumber("Closing iterations", nClosing);
	Dialog.addCheckbox("Detect in soma?", bSomaDetect);
	Dialog.addMessage("Saving options");
	Dialog.addCheckbox("Save in image folder?", bFolder);
	Dialog.addChoice("Save as", newArray("cs and cell ID", "fullname"), saveName);
	Dialog.addCheckbox("Visual options", false);
	Dialog.show;
	mC = Dialog.getNumber();
	sC = Dialog.getNumber();
	xC = Dialog.getNumber();
	minSoma = Dialog.getNumber();
	maxSoma = Dialog.getNumber();
	oneSoma = Dialog.getCheckbox();
	pxSat = Dialog.getNumber();
	gbS = Dialog.getNumber();
	lineW = Dialog.getNumber();
	minLength = Dialog.getNumber();
	highC = Dialog.getNumber();
	lowC = Dialog.getNumber();
	wth = Dialog.getNumber();
	autoT = Dialog.getNumber();
	synS = Dialog.getNumber();
	nThr = Dialog.getNumber();
	neuritePad = Dialog.getNumber();
	nClosing = Dialog.getNumber();
	bSomaDetect = Dialog.getCheckbox();
	bFolder = Dialog.getCheckbox();
	saveName = Dialog.getChoice();
	bVisual = Dialog.getCheckbox();
	if (bVisual) {
		allColor = newArray("white", "orange", "yellow", "cyan", "magenta", "red", "blue", "green");
		Dialog.create("SynJ - visual options");
		Dialog.addChoice("Morphology color", allColor, morRoiColor);
		Dialog.addChoice("Keep synapses color", allColor, keepSynColor);
		Dialog.addChoice("Delete synapses color", allColor, delSynColor);
		Dialog.show();
		morRoiColor = Dialog.getChoice();
		keepSynColor = Dialog.getChoice();
		delSynColor = Dialog.getChoice();
	}
	writeOptions();
}

function autoDetectSoma() {
	imgID = getImageID();
	Stack.setDisplayMode("color");
	Stack.setChannel(mC);
	setBatchMode(true);
	run("Duplicate...", " ");
	dupID = getImageID();
	run("Duplicate...", " ");
	gssID = getImageID();
	run("Minimum...", "radius="+maxSoma);
	run("Maximum...", "radius="+maxSoma);
	run("Gaussian Blur...", "sigma="+maxSoma);
	resetMinAndMax();
	imageCalculator("Subtract", dupID, gssID);
	selectImage(gssID); close();
	run("Minimum...", "radius="+minSoma);
	run("Maximum...", "radius="+minSoma);
	run("Gaussian Blur...", "sigma="+minSoma);
	resetMinAndMax();
	setAutoThreshold("MaxEntropy dark no-reset");
	if (oneSoma) {
		run("Analyze Particles...", "size="+pow(minSoma,2)+"-Infinity pixel add");
		close();
		// filter the biggest most Roi
		nRoi = roiManager("count");
		roiArea = newArray(nRoi);
		for (r = 0; r < nRoi; r++) {
			roiManager("select", r);
			getStatistics(roiArea[r], mean, min, max, std, histogram);
		}
		ranks = Array.rankPositions(roiArea);
		roiManager("Select", ranks[nRoi-1]);
		roiManager("Reset");
		newImage("Untitled", "8-bit black", 1024, 1024, 1);
		mskID = getImageID();
		run("Restore Selection");
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		imageCalculator("AND", mskID, imgID);
		roiManager("Show None");
		run("Gaussian Blur...", "sigma=2");
		setAutoThreshold("MinError dark no-reset");
		run("Convert to Mask");
		run("Erode");
		run("Create Selection");
		selectImage(imgID);
		run("Restore Selection");
		selectImage(mskID); close();
	} else {
		run("Create Selection");
		close();
		run("Restore Selection");
	}
	roiManager("add");
	roiManager("select", 0);
	roiManager("rename", "Soma");
	setBatchMode("Exit and Display");
	call("ij.gui.Toolbar.setBrushSize", minSoma);
}

function editSoma() {
	setTool("brush");
	roiManager("Show None");
	roiManager("select", 0);
	setSlice(mC);
	waitForUser("Manually edit soma");
	roiManager("Update");
}

function editNeurites() {
	setTool("brush");
	roiManager("Show None");
	roiManager("select", 1);
	setSlice(mC);
	waitForUser("Manually edit neurites");
	roiManager("Update");
}

function autoDetectNeurite(){
	Stack.setDisplayMode("color");
	Stack.setChannel(mC);
	// first calculate the threshold for the ridge detection
	rdS = ((lineW)/(2*sqrt(3)))+0.5;
	w2 = (lineW)/2;
	s2 = pow(rdS,2);
	s3 = pow(rdS,3);
	w22 = pow(w2,2);
	exp1 = w22/(2*s2);
	sqrt2pi = sqrt(2*PI);
	denominator = sqrt2pi * s3;
	highN = 2 * highC * w2;
	lowN = 2 * lowC * w2;
	exp2 = exp(-exp1);
	highT = 0.17 * (highN/denominator) * exp2;
	lowT = 0.17 * (lowN/denominator) * exp2;
	setBatchMode(true);
	roiManager("deselect");
	imgID = getImageID();
	run("Duplicate...", "duplicate channels="+mC);
	neuritesID = getImageID();
	//try to clean the background with gaussian smoothing
	run("Duplicate...", " ");
	bkgID = getImageID();
	run("Gaussian Blur...", "sigma="+autoT);
	imageCalculator("Subtract", neuritesID,bkgID);
	selectImage(bkgID); close();
	getMinAndMax(fmin, fmax);
	setMinAndMax(fmin, (fmax/3));
	run("8-bit");
	run("Gaussian Blur...", "sigma="+gbS);
	run("Ridge Detection", "line_width="+lineW+" high_contrast="+highC+" low_contrast="+lowC+" estimate_width extend_line make_binary method_for_overlap_resolution=SLOPE sigma="+rdS+" lower_threshold="+lowT+" upper_threshold="+highT+" minimum_line_length="+minLength+" maximum=0");
	maskID = getImageID();
	//selectImage(imgID);
	//run("Remove Overlay");
	//selectImage(maskID);
	roiManager("select", 0);
	setForegroundColor(255, 255, 255);
	roiManager("fill");
	run("Invert LUT");
	for (c = 0; c < nClosing; c++) {
		run("Close-");
	}
	run("Create Selection");
	selectImage(imgID);
	run("Restore Selection");
	roiManager("add");
	roiManager("select", 1);
	roiManager("rename", "Neurites");
	selectImage(neuritesID); close();
	selectImage(maskID); close();
	setBatchMode("Exit and Display");
	call("ij.gui.Toolbar.setBrushSize", lineW); 
}

function autoDetectSynapses(){
	Stack.setDisplayMode("color");
	Stack.setChannel(sC);
	getPixelSize(unit, pixelWidth, pixelHeight);
	nIter = round(neuritePad/pixelWidth);
	// first calculate the threshold for the ridge detection
	run("Select None");
	imgID = getImageID();
	setBatchMode("hide");
	run("Duplicate...", "duplicate channels="+sC);
	synapsesID = getImageID();
	resetMinAndMax();
	run("Duplicate...", " ");
	bkgID = getImageID();
	run("Gaussian Blur...", "sigma="+autoT);
	imageCalculator("Subtract", synapsesID,bkgID);
	selectImage(bkgID); close();
	run("Morphological Filters", "operation=[White Top Hat] element=Square radius="+wth);
	maskID = getImageID();
	run("Gaussian Blur...", "sigma="+synS);
	roiManager("Select", 1);
	//runMacro(getDirectory("macros") + "autoLUT.ijm"); 
	getMinAndMax(fmin, fmax);
	setMinAndMax(fmin, (fmax/3));
	run("8-bit");
	run("Auto Local Threshold", "method=MidGrey radius="+autoT+" parameter_1=0 parameter_2=0 white");
	getDimensions(width, height, channels, slices, frames);
	newImage("tempMask", "8-bit black", width, height, 1);
	if (bSomaDetect) {
		roiManager("select", 0);
		run("Fill");
	}
	roiManager("select", 1);
	run("Fill");
	run("Select None");
	for (i = 1; i <= nIter; i++)
		run("Dilate");
	setAutoThreshold("Default dark no-reset");
	run("Create Selection");
	selectImage(maskID);
	run("Select None");
	selectWindow("tempMask");
	selectImage(maskID);
	run("Restore Selection");
	selectWindow("tempMask"); close();
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	run("Select None");
	run("Dilate");
	run("Watershed");
	run("Watershed");
	run("Analyze Particles...", "size="+autoT+"-Infinity pixel add");
	selectImage(synapsesID); close();
	selectImage(maskID); close();
	// apply a threshold to filter Roi
	roiManager("select", 1);
	Stack.setChannel(sC);
	getStatistics(nArea, nMean, nMin, nMax, nStd, nHistogram);
	if (nThr != 0) {
		threshold = nMean + nThr * nStd;
		for (r = 2; r < roiManager("count"); r++) {
			roiManager("select", r);
			roiManager("rename", "Synapse_"+r+1);
			Stack.setChannel(sC);
			getStatistics(rArea, rMean, rMin, rMax, rStd, rHistogram);
			if (rMean <= threshold) {
				Roi.setStrokeColor(delSynColor);
			} else {
				Roi.setStrokeColor(keepSynColor);
				}
		}
	} else {
		for (r = 2; r < roiManager("count"); r++) {
			roiManager("select", r);
			Roi.setStrokeColor(keepSynColor);
		}
	}
	setBatchMode("Exit and display");
}

function deleteWrongSynapses() {
	// when the synapses are selected than delete them
	setBatchMode("hide");
	nRoi = roiManager("count");
	for (r = nRoi-1; r >= 0; r--) {
		showProgress(nRoi-r,nRoi);
		showStatus("Cleaning Roi Manager " +r+"/"+nRoi);
		roiManager("Select", r);
		bDel = matches(Roi.getStrokeColor(), delSynColor);
		if(bDel){
			roiManager("Delete");
		}
	}
	/*
	r = 0;
	while(r<nRoi){
		showProgress(-(r+1),nRoi);
		showStatus("Cleaning Roi Manager " +r+"/"+nRoi);
		roiManager("Select", r);
		r1 = 0;
		bDel = matches(Roi.getStrokeColor(), delSynColor);
		if(bDel){
			roiManager("Delete");
			nRoi = roiManager("count");
		}else{
			r++;
		}
	}
	*/
	setSlice(sC);
	run("Remove Overlay");
	setBatchMode("Exit and display");
}

function deleteAllSynapses(){
	setBatchMode(true);
	imgID = getImageID();
	getDimensions(width, height, channels, slices, frames);
	newImage("Soma", "8-bit", width, height, 1);
	roiManager("Select", 0);
	newImage("Neurite", "8-bit", width, height, 1);
	roiManager("Select", 1);
	selectImage(imgID);
	roiManager("Reset");
	selectWindow("Soma");
	roiManager("add"); roiManager("Select", 0); roiManager("rename", "Soma");
	close();
	selectWindow("Neurite");
	roiManager("add"); roiManager("Select", 1); roiManager("rename", "Neurite");
	close();
	//nRoi = roiManager("count");
	//r = 0;
	//while (r<nRoi) {
	//	showProgress(-(r+1),nRoi);
	//	showStatus("Cleaning Roi Manager " +r+"/"+nRoi);
	//	if (roiManager("count") > 2) {
	//		roiManager("select", r);
	//		if (startsWith(Roi.getName,"Synapse_")) {
	//			roiManager("delete");
	//		} else {
	//			r++;
	//		}
	//	}
	//}
}

function saveIntensity(saveDir,nameToSave) {
	// first calculate the skeleton lenght since it's rather tedius
	getDimensions(width, height, channels, slices, frames);
	setBatchMode(true);
	newImage("tempSkeleton", "8-bit black", width, height, 1);
	roiManager("select", 1);
	run("Fill");
	run("Select None");
	setOption("BlackBackground", true);
	run("Dilate");
	run("Options...", "iterations=1 count=4 black do=Nothing"); // to smooth the neurite mask
	run("Erode");
	run("Skeletonize");
	run("Analyze Skeleton (2D/3D)", "prune=[shortest branch]");
	run("Options...", "iterations=1 count=1 black do=Nothing");
	totLength = 0;
	endPoint = 0;
	for (s = 0; s < nResults; s++) {
		totLength = totLength + getResult("# Branches", s) * getResult("Average Branch Length", s);
		endPoint = endPoint + getResult("# End-point voxels", s);
	}
	selectWindow("Results"); run("Close");
	selectWindow("tempSkeleton"); run("Close");
	if (isOpen("Tagged skeleton")) {
		selectWindow("Tagged skeleton"); run("Close");
	}
	//setBatchMode(false);
	// get the value for each roi, not only synapses, but also soma and neurites
	setBatchMode("Hide");
	nRoi = roiManager("count");
	run("Clear Results");
	res = 0;
	for (r = 0; r < nRoi; r++) {
		showProgress(r/nRoi);
		showStatus("Saving Intensity");
		roiManager("select", r);
		if (matches(Roi.getStrokeColor(), delSynColor)) {
			
		} else {
			Stack.setChannel(mC);
			getStatistics(roiArea, roiMean, roiMin, roiMax, roiStd, roiHistogram);
			if (r == 1) {
				setResult("Region", res, Roi.getName+" ("+endPoint+")");
				setResult("Area", res, totLength);
			} else {
				setResult("Region", res, Roi.getName);
				setResult("Area", res, roiArea);
			}
			setResult("Morphology", res, roiMean);
			Stack.setChannel(sC);
			getStatistics(roiArea, roiMean, roiMin, roiMax, roiStd, roiHistogram);
			setResult("Synapses", res, roiMean);
			Stack.setChannel(xC);
			getStatistics(roiArea, roiMean, roiMin, roiMax, roiStd, roiHistogram);
			setResult("Other", res, roiMean);
			res++;
		}
	}
	updateResults();
	saveAs("Results_"+nameToSave, saveDir + "/Results_" + nameToSave + ".csv");
	selectWindow("Results"); run("Close");
	setBatchMode("Exit and display");
}

function SynJ_Sholl(regSholl, sRadius) {
	// Sholl intensity analysis
	imgID = getImageID();
	Stack.setDisplayMode("color");
	// Define center
	roiManager("Select", 0);
	Roi.getBounds(xSoma, ySoma, wSoma, hSoma);
	xCenter = xSoma + wSoma / 2;
	yCenter = ySoma + hSoma / 2;
	run("Select None");
	// get Radius (converted in pixel)
	getPixelSize(unit, wPx, hPx);
	sRadius = round(sRadius / wPx);
	stRadius = maxOf(wSoma/2, hSoma/2);
	roiManager("Select", 1);
	Roi.getContainedPoints(xNeu, yNeu);
	endRadius1 = sqrt(pow(xNeu[0] - xCenter, 2) + pow(yNeu[0] - yCenter, 2));
	endRadius2 = sqrt(pow(xNeu[lengthOf(xNeu)-1] - xCenter, 2) + pow(yNeu[lengthOf(yNeu)-1] - yCenter, 2));
	endRadius = maxOf(endRadius1, endRadius2);
	nRounds = round((endRadius - stRadius) / sRadius) + 1;
	// create a new image for the actual region to measure
	setBatchMode(true);
	newImage("Dummy_Neurite", "8-bit black", getWidth(), getHeight(), 1);
	neuriteID = getImageID();
	roiManager("Select", 1);
	roiManager("fill");
	newImage("Dummy_Synapses", "8-bit black", getWidth(), getHeight(), 1);
	nRoi = roiManager("count");
	for (r = 0; r < nRoi; r++) {
		roiManager("Select", r);
		if (matches(Roi.getStrokeColor(), keepSynColor)) {
			roiManager("Fill");
		}
	}
	synapsesID = getImageID();
	newImage("Dummy", "8-bit black", getWidth(), getHeight(), 1);
	dumID = getImageID();
	// start makind the circles, first get the soma intensity
	selectImage(imgID);
	setResult("Distance", 0, "Soma");
	roiManager("Select", 0);
	Stack.setChannel(mC);
	getStatistics(sArea, mMean, sMin, sMax, sStd, sHistogram);
	setResult("Morphology(N)", 0, mMean);
	Stack.setChannel(sC);
	getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
	setResult("Synapses(N)", 0, sMean);
	Stack.setChannel(xC);
	getStatistics(sArea, xMean, sMin, sMax, sStd, sHistogram);
	setResult("Other(N)", 0, xMean);
	for (r=0; r<=nRounds; r++) {
		selectImage(dumID);
		run("Select All");
		run("Clear");
		thisR = stRadius + r * sRadius;
		setResult("Distance", r+1, thisR*wPx);
		makeOval(xCenter-thisR, yCenter-thisR, thisR*2, thisR*2);
		run("Make Inverse");
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		makeOval(xCenter-(thisR+sRadius), yCenter-(thisR+sRadius), (thisR+sRadius)*2, (thisR+sRadius)*2);
		run("Clear Outside");
		run("Select None");
		if (matches(regSholl, "Neurite")) {
			imageCalculator("AND create", dumID, neuriteID);
			tempNeuID = getImageID();
			run("Make Binary");
			run("Create Selection");
			if (selectionType()!=-1) {
				selectImage(imgID);
				run("Restore Selection");
				Stack.setChannel(mC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Morphology(N)", r+1, sMean);
				Stack.setChannel(sC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Synapses(N)", r+1, sMean);
				Stack.setChannel(xC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Other(N)", r+1, sMean);
			} else {
				setResult("Morphology(N)", r+1, "NaN");
				setResult("Synapses(N)", r+1, "NaN");
				setResult("Other(N)", r+1, "NaN");
			}
			selectImage(tempNeuID); close();
		} else if (matches(regSholl, "Synapses")) {
			imageCalculator("AND create", dumID, synapsesID);
			tempSynID = getImageID();
			run("Make Binary");
			run("Create Selection");
			if (selectionType()!=-1) {
				selectImage(imgID);
				run("Restore Selection");
				Stack.setChannel(mC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Morphology(S)", r+1, sMean);
				Stack.setChannel(sC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Synapses(S)", r+1, sMean);
				Stack.setChannel(xC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Other(S)", r+1, sMean);
			} else {
				setResult("Morphology(S)", r+1, "NaN");
				setResult("Synapses(S)", r+1, "NaN");
				setResult("Other(S)", r+1, "NaN");
			}
			selectImage(tempSynID); close();
		} else {
			imageCalculator("AND create", dumID, neuriteID);
			tempNeuID = getImageID();
			run("Make Binary");
			run("Create Selection");
			if (selectionType()!=-1) {
				selectImage(imgID);
				run("Restore Selection");
				Stack.setChannel(mC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Morphology(N)", r+1, sMean);
				Stack.setChannel(sC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Synapses(N)", r+1, sMean);
				Stack.setChannel(xC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Other(N)", r+1, sMean);
			} else {
				setResult("Morphology(N)", r+1, "NaN");
				setResult("Synapses(N)", r+1, "NaN");
				setResult("Other(N)", r+1, "NaN");
			}
			// synapses
			selectImage(tempNeuID); close();
			imageCalculator("AND create", dumID, synapsesID);
			tempSynID = getImageID();
			run("Make Binary");
			run("Create Selection");
			if (selectionType()!=-1) {
				selectImage(imgID);
				run("Restore Selection");
				Stack.setChannel(mC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Morphology(S)", r+1, sMean);
				Stack.setChannel(sC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Synapses(S)", r+1, sMean);
				Stack.setChannel(xC);
				getStatistics(sArea, sMean, sMin, sMax, sStd, sHistogram);
				setResult("Other(S)", r+1, sMean);
				selectImage(tempSynID); close();
			} else {
				setResult("Morphology(S)", r+1, "NaN");
				setResult("Synapses(S)", r+1, "NaN");
				setResult("Other(S)", r+1, "NaN");
			}
		}
	}
	selectImage(dumID); close();
	selectImage(neuriteID); close();
	selectImage(synapsesID); close();
	updateResults();
}

function getSaveName() {
		// get image folder
	if (bFolder){
  		saveDir = getDirectory("Image");
  	} else {
		saveDir = getDirectory("Choose a Directory");
  	}
  	
  	// get the name
	imgName = getTitle();
	fileName = split(imgName, ".");
	imgName = fileName[0];
	
	// get how to save
	if (saveName == "cs and cell ID"){
		namePart = split(imgName, "_");
		nameToSave = namePart[3]+"_"+namePart[4];
	} else {
		nameToSave = imgName;
	}
}

function adjustTable() {
	imgDir = getDirectory("Select images folder");
	imgFile = getFileList(imgDir);
	nImg = imgFile.length;
	for (i = 0; i < nImg; i++) { // loop throught the images
		if (endsWith(imgFile[i], ".tif")) {
			for (r = 0; r < nImg; r++) { // loop to find the right RoiSet
				bRoi = startsWith(imgFile[r], "RoiSet");
				bImg = indexOf(imgFile[r], substring(imgFile[i],0,lengthOf(imgFile[i])-4)) > 0;
				if (bImg && bRoi) {
					setBatchMode(true);
					open(imgDir + imgFile[i]);
					orTitle = getTitle();
					orTitle = replace(orTitle, ".tif", "");
					getDimensions(width, height, channels, slices, frames);
					imgID = getImageID();
					roiManager("Open",  imgDir + imgFile[r]);
					newImage("tempSkeleton", "8-bit black", width, height, 1);
					roiManager("select", 1);
					run("Fill");
					run("Select None");
					setOption("BlackBackground", true);
					run("Dilate");
					run("Options...", "iterations=1 count=4 black do=Nothing"); // to smooth the neurite mask
					run("Erode");
					run("Skeletonize");
					run("Analyze Skeleton (2D/3D)", "prune=[shortest branch]");
					run("Options...", "iterations=1 count=1 black do=Nothing");
					totLength = 0;
					endPoint = 0;
					for (s = 0; s < nResults; s++) {
						totLength = totLength + getResult("# Branches", s) * getResult("Average Branch Length", s);
						endPoint = endPoint + getResult("# End-point voxels", s);
					}
					selectWindow("Results"); run("Close");
					selectWindow("tempSkeleton"); run("Close");
					if (isOpen("Tagged skeleton")) {
						selectWindow("Tagged skeleton"); run("Close");
					}
					//setBatchMode(false);
					open(imgDir + "/Results_" + orTitle + ".csv");
					Table.rename("Results_"+orTitle+".csv", "Results");
					setResult("Region", 1, "Neurites ("+endPoint+")");
					updateResults();
					saveAs("Results_"+orTitle, imgDir + "/Results_" + orTitle + ".csv");
					selectWindow("Results"); run("Close");
					roiManager("reset");
					selectImage(imgID); close();
				}
			}
		}
	}
}