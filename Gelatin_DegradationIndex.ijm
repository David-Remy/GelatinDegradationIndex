/*  Works on ImageJ 1.53c / Windows
 *  Created in March 2021 by David Remy
 *  Curated by Anne-Sophie Macé 
 *  This macro quantifies gelatin degradation per field and normalize to number of nuclei in the field.
 *  Batch quantification included: input should be a folder containing all acquisitions
*/

 // Clears everything
roiManager("reset");
print("\\Clear");
run("Clear Results");
run("Close All");
run("Set Measurements...", "area centroid stack redirect=None decimal=3");
setOption("ExpandableArrays", true); // In ImageJ 1.53g and later, arrays automatically expand in size as needed. This option is used in case of early versions of ImageJ
if( isOpen("Summary") ){
	selectWindow("Summary");
	run("Close");
}

/////////////////////////////////////////////////////////////
////// begining of parameters customozible by the user //////
/////////////////////////////////////////////////////////////
// extensions of the file to study [tested on ND only]
ext_file = ".nd";
// minimum area in pixels of the nuclei to detect
nuclei_min_area = 1000;
// size of the radius in pixels for the unsharp mask of gelatin 
unsharp_radius = 5;
// Mask size in pixels for the unsharp mask of gelatin 
unsharp_weight = 0.60;
// size of the radius size in pixels for the mean mask of gelatin
mean_radius = 0.80;
// minimum area in pixels for Analyze particle in the gelatin degradation detection
analyze_particle_min_gel = 10;
/////////////////////////////////////////////////////////////
//////// end of parameters customozible by the user /////////
/////////////////////////////////////////////////////////////
 

// Select input directory and create an Analysis subfolder
dir_input=getDirectory("Select input directory");
if( !File.exists(dir_input+"Gelatin_Degradation_Analysis") ) { // creates directory if does not exist
	File.makeDirectory(dir_input+"Gelatin_Degradation_Analysis");
}
dir_output=dir_input+"Gelatin_Degradation_Analysis"+File.separator; 

// Dialog box for channel number/name assignment 
channels_array=newArray(1, 2, 3, 4);
Dialog.create("Channel representation");
Dialog.addMessage("1: DAPI ; 2:GFP ; 3: Cy3; 4: Cy5"); 
Dialog.addChoice("Nucleus channel", channels_array, "1.0");
Dialog.addChoice("Gelatin channel", channels_array, "2.0");

Dialog.addCheckbox("If z-stack, Max Projection on the entire z-stack ?" , false);
Dialog.show();

Nuc_channel=Dialog.getChoice();
Degra_channel=Dialog.getChoice();
proj=Dialog.getCheckbox();

//list every files names within the directory
Filelist=getFileList(dir_input);
Array.sort(Filelist);

//Initializes arrays to compile all names of images, number of nuclei per image, area of degradation per image and area/numb of nuclei
Name=newArray();
NumNuclei=newArray();
AreaDegradation=newArray();
normArea=newArray();

nd_file=0;
file_treated = 0;
for (i_file=0; i_file<lengthOf(Filelist); i_file++) {
	if(indexOf(Filelist[i_file], ext_file)>0){ // treats each ND file
		shortTitle = substring(Filelist[i_file],0,lastIndexOf(Filelist[i_file],"."));
		
		run("Bio-Formats Importer", "open=["+dir_input+Filelist[i_file]+"] autoscale color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		tit_img = getTitle();
		
		// convert scale in pixels
		run("Set Scale...", "distance=1");
		getDimensions(width, height, channels, slices, frames);
		
		// check that the specified channels exist
		if( channels < maxOf(Degra_channel,Nuc_channel) )
			print("Image "+Filelist[i_file]+" not treated: at least one of the specified channel do not exist");
		
		else{
			print("Treating image "+Filelist[i_file]);
			
			if( !File.exists(dir_output+"Results_"+shortTitle+".xls") ){ // the image was not treated
				// duplicates the channels of interest
				selectWindow(tit_img);
				run("Duplicate...", "title=gelatin duplicate channels="+Degra_channel);
				selectWindow(tit_img);
				run("Duplicate...", "title=nucleus duplicate channels="+Nuc_channel);
				
				//counts the number of nuclei in the field on the z-max projection
				selectWindow("nucleus");	
				if (nSlices >1 ) { // If the image is a stack: maximum z-projection
					run("Z Project...", "projection=[Max Intensity]");
					close("nucleus");
					selectWindow("MAX_nucleus");
					rename("nucleus");
				}
				run("Gaussian Blur...", "sigma=2");
				setAutoThreshold("Default dark");
				run("Analyze Particles...", "size="+nuclei_min_area+"-Infinity exclude display"); //The "exclude" command excludes nuclei that are on the borders of the image
				
				// Nuclei = number of found particules; a dialog box asks the user to confirm or modify this number
				nuclei=nResults;
				NumNuclei[file_treated]=getNumber("Nuclei detected", nuclei);
				
				selectWindow("gelatin");	
				first_plan = 1;
				last_plan = nSlices;
				
				// image is a stack and the user asked to choose the plans of analysis
				if( nSlices > 1 && !proj) {
					// asks the user the z plans on which the projection should be performed
					Dialog.createNonBlocking("Plans choice");
					Dialog.addMessage("Choose the plan to analyze, between 1 and "+slices+" (put the same value if you want only one image)."); 
					Dialog.addMessage("If several plans, a projection will be performed.");
					Dialog.addNumber("First plan (above 1)", first_plan);
					Dialog.addNumber("Last Plan (below "+slices+")", last_plan);
					Dialog.show();
				
					first_plan = Dialog.getNumber();
					last_plan = Dialog.getNumber();
					
					if( first_plan < 1 || first_plan > slices || last_plan < 1 || last_plan > slices)
						print("You choose an uncorrect number for first or last plan, the maximum projection is computed on all slices");
			
					// apply maximum z-projection on the slices specified by the user
					selectWindow("gelatin");
					run("Z Project...", "start="+first_plan+" stop="+last_plan+" projection=[Max Intensity]");
					close("gelatin");
				}
				// Max Projection on all the stack if the user chose that option 
				if( nSlices > 1 && proj) 
					run("Z Project...", "start="+first_plan+" stop="+last_plan+" projection=[Max Intensity]");	
				else // to have the same name
					rename("MAX_gelatin");
				
				selectWindow("MAX_gelatin");
				
				// pre-processing - user defined threshold - analyze particle
				run("Unsharp Mask...", "radius="+unsharp_radius+" mask="+unsharp_weight);
				run("Mean...", "radius="+mean_radius);
				run("Threshold...");
				setThreshold(0, 1358, "raw");
				waitForUser("Please select the appropriate thresold. You can draw a delimiting region if required. Then, press OK !");
				run("Analyze Particles...", "size="+analyze_particle_min_gel+"-Infinity clear summarize ");
				
				// saves updated summaries and thresholds 
				close("*");
				selectWindow("Summary");
				IJ.renameResults("Summary","Results");
				Name[file_treated]=shortTitle;
				AreaDegradation[file_treated]=getResult("Total Area", 0);
				normArea[file_treated]=AreaDegradation[file_treated]/NumNuclei[file_treated];
				setResult("NucleiNumber",0,NumNuclei[file_treated]); // to be able to read again results if the macro crashes/the user had to stop
				saveAs("Results", dir_output+"Results_"+shortTitle+".xls");
				file_treated++;
				close("Results");
			}
			else{ // the image was aready treated: we load the results
				print("Results file existed: loaded");
				run("Results... ","open=["+dir_output+"Results_"+shortTitle+".xls]");
				Name[file_treated]=shortTitle;
				AreaDegradation[file_treated] = getResult("Total Area", 0);
				NumNuclei[file_treated] = getResult("NucleiNumber", 0);
				normArea[file_treated] = AreaDegradation[file_treated]/NumNuclei[file_treated];
				file_treated++;
				close("Results");
			}
		}
	}	
}

// final result table for all images of the folder
run("Clear Results");
for (i_results = 0; i_results < lengthOf(NumNuclei); i_results++) {
	setResult("Name", i_results, Name[i_results]);
	setResult("Nuclei", i_results, NumNuclei[i_results]);
	setResult("Total Area", i_results, AreaDegradation[i_results]);
	setResult("Area/Nuclei", i_results, normArea[i_results]);
}	
saveAs("Results", dir_output+File.getName(dir_input)+"_GelatinNormalizedArea.xls");
run("Close All");