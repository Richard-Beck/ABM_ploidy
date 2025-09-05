package onLatticeCA_jar;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.UIGrid;
import HAL.Tools.FileIO;
import HAL.Util;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.sql.SQLOutput;
import java.util.Arrays;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.io.*;

import static HAL.Util.*;


public class ExportData{
    /**
     * For exporting ABM data for further analysis
     */
    File outputFolder,imageFolder,karyotypeFolder,oxygenFolder;
    FileIO oxygenSummary,summary;
    public int karyotypeDrawValue;
    public String karyotypeWriteValue;
    public OnLatticeGrid G;
   // public OnLatticeGrid G;
    Params params=new Params();
    helperMethods helper=new helperMethods();
    Resources resources=new Resources();


    public ExportData(String outputPath, String path2Params) throws IOException {
        // we make a new folder each time we run the model.
        // the folder will be a subdirectory of the outputPath argument.
        File outputHome = new File(outputPath);
        File[] folderList = outputHome.listFiles();
        outputPath = outputPath + "/run_" + folderList.length;
        outputFolder = new File(outputPath);
        outputFolder.mkdirs();
        
        //make a copy of the parameter file in the output directory:
        Path paramSource = Paths.get(path2Params);
        Path paramDest = Paths.get(outputPath + "/params.txt");
        Files.copy(paramSource, paramDest, StandardCopyOption.REPLACE_EXISTING);

        // create subfolders to store any images.
        imageFolder = new File(outputPath+"/images");
        imageFolder.mkdirs();

        karyotypeFolder = new File(outputPath+"/karyotypes");
        karyotypeFolder.mkdirs();

        oxygenFolder = new File(outputPath+"/oxygen");
        oxygenFolder.mkdirs();
        // create any files which we intend to write to throughout the simulation
        oxygenSummary = new FileIO(outputPath + "/oxygen.csv", "w");
        oxygenSummary.Write("time, O2, maxDelta \n");

        summary = new FileIO(outputPath + "/summary.csv", "w");
        summary.Write("time,nNormal,Ncancer \n");
    }
    public ExportData(){}

    public void saveImage(UIGrid cells, UIGrid vessels, int tstep){
        UIGrid vis = new UIGrid(cells.xDim+vessels.xDim, cells.yDim, 2);
        for(int y=0; y<cells.yDim; y++){
            for(int x=0; x<cells.xDim; x++){
                vis.SetPix(x,y,cells.GetPix(x,y));
                vis.SetPix(x+cells.xDim,y,vessels.GetPix(x,y));
            }
        }
        String outputPath = imageFolder +"/"+tstep+".png";
        vis.ToPNG(outputPath);
    }


    public void writeOxygenSummary(double time, double O2, double maxDelta){
        oxygenSummary.Write(time+"," + O2 +"," + maxDelta + "\n");
    }
    public void writeSummary(double time, int nNormal, int Ncancer){
        summary.Write(time+"," + nNormal+ "," + Ncancer + "\n");
    }
    public void close(){

        summary.Close();
        oxygenSummary.Close();
    }
    public void saveKaryotypeLocationData(int currentTime,OnLatticeGrid G){
        FileIO karyotype_Location = new FileIO(karyotypeFolder+"/"+currentTime +".csv","w");
        for (Cell c : G) {
            int k = c.Isq();
            karyotype_Location.Write(G.ItoX(k) +","+ G.ItoY(k));
            for(int m=0;m<22;m++)  karyotype_Location.Write(","+c.karyotype[m]);
            karyotype_Location.Write("\n");
        }
        karyotype_Location.Close();
    }
    public void saveOxygenField(int currentTime, PDEGrid2D o2){
        FileIO oxygen_Location = new FileIO(oxygenFolder+"/"+currentTime +".csv","w");
        for(int y=0;y<o2.yDim;y++){
            oxygen_Location.Write(o2.Get(0,y)+"");
            for(int x=1;x<o2.xDim;x++){
                oxygen_Location.Write(","+o2.Get(x,y));
            }
            oxygen_Location.Write("\n");
        }
        oxygen_Location.Close();
    }

    // Call this once per tick, handles images + heavy dumps by cadence
    public void maybeSnapshot(int tIdx, UIGrid vis, Resources res, OnLatticeGrid model, double dt){
        boolean wantImg = (tIdx == 1) || (Params.imageFrequency > 0 && (tIdx % (int)(Params.imageFrequency/dt) == 0));
        if (wantImg) saveImage(vis, res.currV, tIdx);

        boolean wantData = Params.writeKaryotypeData &&
                        ((tIdx == 1) || (Params.writeFrequency > 0 && (tIdx % (int)(Params.writeFrequency/dt) == 0)));
        if (wantData) {
            saveKaryotypeLocationData(tIdx, model);
            saveOxygenField(tIdx, res.pdegrid2d);
        }
    }

    // One place to log the per-tick summary + oxygen summary
    public void writeTick(int tIdx, int normal, int cancer, double simTimeAdded, double meanO2, double maxDelta){
        writeSummary((double)tIdx, normal, cancer);
        writeOxygenSummary(simTimeAdded, meanO2, maxDelta);
    }
}