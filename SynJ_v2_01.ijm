// SynJ - fluoresce synapse analysis in ImageJ

/*

Start developing 2018.10.21
Start release 2019.11.30

Modify
	20.01.21 - Start creating v 2.0 with separate functions
	20.07.14 - Implemented detection in two channels and colocalization (See separate functions), some cleaning and bug fixing  under the hood
*/


var majVer = 2;
var minVer = 01;
var about = "Developed by Alessandro Moro<br>"
			+ "<i>Department of Functional Genomics</i> (FGA)<br>"
			+ "<i>Centre of neuroscience and cognitive research</i> (CNCR)<br>"
			+ "<i>Vrij Universiteit</i> (VU) Amsterdam.<br>"
			+ "<i>email: a.moro@vu.nl</i><br><br><br>";

var SynJ_dir = getDirectory("imagej") + "macros//toolsets//SynJ";
// even before starting check if it's the first time it's run
function firstCheck() {
	bFirst = call("ij.Prefs.get", "SynJ.firstRun", true);
	if(bFirst == 1){
		runMacro(SynJ_dir+"//SynJ_Options.ijm", "firstTime"); 
		call("ij.Prefs.set", "SynJ.firstRun", false);
	}
}

/////////////////////////////////////
/////////CREATE THE TOOLSET/////////
////////////////////////////////////

// leave one empty slot
macro "Unused Tool -1-" {}

// Soma Analysis -> Automatic detection of soma? Define soma mask
var sCmds1 = newMenu("Detect soma Menu Tool", newArray("Auto detection", "Define soma", "Edit soma", "-", "Get from SynD", "Batch analysis", "-", "Convert *.nd2"));
macro "Detect soma Menu Tool - Ca00 T0b10S T6b10o Tcb10m"{
	firstCheck();
	cmd1 = getArgument();
	if (cmd1 == "Auto detection") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_SomaDetection.ijm", "autoDetection");
	} else if (cmd1 == "Define soma") {
		//setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_SomaDetection.ijm", "defineSoma");
	} else if (cmd1 == "Edit Soma") {
		runMacro(SynJ_dir+"//SynJ_SomaDetection.ijm", "editSoma");
	} else if (cmd1 == "Get from SynD") {
		runMacro(SynJ_dir+"//SynJ_BatchProcessing.ijm", "getSynD");
	} else if (cmd1 == "Batch analysis") {
		runMacro(SynJ_dir+"//SynJ_BatchProcessing.ijm", "batchAnalysis");
	} else {
		runMacro(SynJ_dir+"//SynJ_BatchProcessing.ijm", "convertND2");
	}
}

// Neurite Analysis -> Automatic detection of neurites, Define neurites mask
var sCmds2 = newMenu("Detect neurites Menu Tool", newArray("Auto detection", "Edit neurites"));
macro "Detect neurites Menu Tool - C0a0 T0b10N T7b10e Tcb10u"{
	firstCheck();
	cmd2 = getArgument();
	if (cmd2 == "Auto detection") {
		runMacro(SynJ_dir+"//SynJ_NeuriteDetection.ijm", "autoDetection");
		run("Remove Overlay");
	} else if (cmd2 == "Edit neurites") {
		runMacro(SynJ_dir+"//SynJ_NeuriteDetection.ijm", "editNeurites");
	}
}

// Synapse Analysis -> Automatic detection of synapses
var sCmds3 = newMenu("Detect synapses Menu Tool", newArray("Auto detection", "Detect other", "Run colocalization", "-", "Delete synapses", "Delete all synapses"));
macro "Detect synapses Menu Tool - C00a T0b10S T7b10y Tcb10n"{
	firstCheck();
	cmd3 = getArgument();
	if (cmd3 == "Auto detection") {
		runMacro(SynJ_dir+"//SynJ_SynapseDetection.ijm", "autoDetection");
		setBatchMode("Exit and display");
		run("Remove Overlay");
	} else if (cmd3 == "Delete synapses") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_SynapseDetection.ijm", "deleteSynapses");
		setBatchMode("Exit and display");
	} else if (cmd3 == "Detect other") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_SynapseDetection.ijm", "detectPost");
		setBatchMode("Exit and display");
		run("Remove Overlay");
	} else if (cmd3 == "Delete all synapses") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_SynapseDetection.ijm", "deleteAllSynapses");
	} else if (cmd3 == "Run colocalization") {
		runMacro(SynJ_dir+"//SynJ_colocalizationTest2.ijm");
	}
}

// Mark Synapse tool -> toggle between keep and delete color for synapses
macro "Mark syanpses Tool - Cf0f R0055R7255R3a55 Cf11 L33ccL43bcLc33cLb33b" {
	sC = call("ij.Prefs.get", "SynJ.SynapseChannel", true);
	keepSynColor = call("ij.Prefs.get", "SynJ.keepSynapsesColor", true);
	delSynColor = call("ij.Prefs.get", "SynJ.deleteSynapsesColor", true);
	setSlice(sC);
	getCursorLoc(x, y, z, flags);
	idx = roiManager("index");
	if (flags == 48) {
		if (idx > 0) {
			roiManager("select", idx);
			if (selectionContains(x, y)) {
				if (matches(Roi.getStrokeColor(), delSynColor)){
					roiManager("Set Color", keepSynColor);
				} else {
					roiManager("Set Color", delSynColor);
				}
			} else {
				run("Select None");
			}
		} else {
			waitForUser("Try pressing longer");
		}
	}
}

// Save
var sCmds4 = newMenu("Save Menu Tool", newArray("Save all", "Roi Manager", "Intensity", "Adjust table", "-", "Sholl Analysis", "-", "Collapse data"));
macro "Save Menu Tool - C555D11D12D13D14D15D16D17D18D19D1aD1bD1cD21D27D28D29D2aD2bD2cD2dD31D33D35D37D38D39D3aD3bD3cD3dD3eD41D43D45D47D48D4eD51D53D55D57D58D5aD5bD5cD5eD61D63D65D67D68D6aD6bD6cD6eD71D73D75D77D78D7eD81D83D85D87D88D8eD91D93D95D97D98D9eDa1Da3Da5Da7Da8DaeDb1Db7Db8DbeDc1Dc2Dc3Dc4Dc5Dc6Dc7Dc8Dc9DcaDcbDccDcdDce"{
	firstCheck();
	cmd4 = getArgument();
	if (cmd4 == "Save all") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "roiManager");
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "intensity");
		setBatchMode("Exit and display");
	} else if (cmd4 == "Roi Manager") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "roiManager");
		setBatchMode("Exit and display");
	} else if (cmd4 == "Intensity") {
		setBatchMode(true);
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "intensity");
		setBatchMode("Exit and display");
	} else if (cmd4 == "Adjust table") {
		bCont = true;
		while (bCont) {
			runMacro(SynJ_dir+"//SynJ_Save.ijm", "adjustTable");
			bCont = getBoolean("Would you like to adjust another set?");
		}
		return;
	} else if (cmd4 == "Sholl Analysis") {
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "shollAnalysis");
	} else if (cmd4 == "Collapse data") {
		runMacro(SynJ_dir+"//SynJ_Save.ijm", "collapseData");
	}
}

// Options
macro "Options... Action Tool - C77bD3eD4eD5eD6bD6cD6dD7aD89D98Da7Db6Dc6Dd6De4De5D2aD5dDa2Dd5D59D68D69D77D78D86D87D96D1aD1bD1cD29D2bD39D49D4bD4cD4dD58D67D76D85D92D93D94Da1Db1Db2Db4Dc1Dc4Dd4De3D5aD6aD79D88D95D97Da5Da6D19D91D4aD5bDa4Db5D3aD5cDa3Dc5"{
	firstCheck();
	runMacro(SynJ_dir+"//SynJ_Options.ijm", "setOptions");
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