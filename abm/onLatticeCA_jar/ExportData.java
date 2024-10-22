package onLatticeCA_jar;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.UIGrid;
import HAL.Tools.FileIO;
import HAL.Util;

import java.io.IOException;
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
import java.io.BufferedWriter;
import java.io.FileWriter;

import static HAL.Util.*;


public class ExportData{
    /**
     * For exporting ABM data for further analysis
     */
    File outputFolder,imageFolder,karyotypeFolder,oxygenFolder,gifFolder;
    FileIO oxygenSummary,summary;
    GifMaker cells_O2;
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
        gifFolder =new File(outputPath+"/gifs");
        gifFolder.mkdirs();

        karyotypeFolder = new File(outputPath+"/karyotypes");
        karyotypeFolder.mkdirs();

        oxygenFolder = new File(outputPath+"/oxygen");
        oxygenFolder.mkdirs();
        // create any files which we intend to write to throughout the simulation
        oxygenSummary = new FileIO(outputPath + "/oxygen.csv", "w");
        oxygenSummary.Write("time, O2, maxDelta \n");

        cells_O2=new GifMaker(outputPath +"/gifs/cell_oxygen.gif",1,true);


        summary = new FileIO(outputPath + "/summary.csv", "w");
        summary.Write("time,nNormal,nCancer \n");
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
        cells_O2.Close();
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
    public void createGifFiles(UIGrid cells, UIGrid vessels){
        UIGrid vis = new UIGrid(cells.xDim+vessels.xDim, cells.yDim, 2);
        for(int y=0; y<cells.yDim; y++){
            for(int x=0; x<cells.xDim; x++){
                vis.SetPix(x,y,cells.GetPix(x,y));
                vis.SetPix(x+cells.xDim,y,vessels.GetPix(x,y));
            }
        }
        cells_O2.AddFrame(vis);
    }

    public void saveGif(UIGrid frame){

    }

    public void saveData(int currentTime,OnLatticeGrid G){


        ExportData exportData = new ExportData();
       // FileIO karyotypeDist = new FileIO(currentTime+"distribution.csv", "w");
        FileIO karyotype_Location = new FileIO(G.mainDir.concat("/output/"+currentTime+"_karyotype_Location.csv"), "w");

        //FileIO cellLoc=new FileIO(G.mainDir.concat("/output/cellLoc.txt"), "w");

        List<int[][]> k1=new ArrayList<>();


        int[][] outputData= new int[G.xDim][G.yDim];

        //int[][]   cellLocations=new int[G.xDim*G.yDim][2];
        int[][] karyotypeAndLocation = new int[(G.xDim * G.yDim)][24];
        double [][] O2ConcentrationPerGridCell = new double[G.xDim][G.yDim];
        double [][] TotalO2ConcentrationInTheGrid=new double[(int)G.tEnd][2];//



        //int[][] Karyotype1= new int[xDim][yDim];
        for (Cell c : G) {
            int k = c.Isq();
            int i = G.ItoX(k); //Converts an index of grid X and Y coordinate
            int j = G.ItoY(k); //Converts an index of grid X and Y coordinate
            //@TODO UPDATE DATA AQUISITION Process



            //Saving 2D array (x,y,z) where x,y are coordinates  z is karyotype related
            karyotypeAndLocation[k][0]=i;
            karyotypeAndLocation[k][1] = j;
            for(int m=0;m<22;m++){
                karyotypeAndLocation[k][m+2]=c.karyotype[m];
            }



        }

       //Write Karyotype-location files
        /**
        karyotypeDist.WriteDelimit(outputData[0], ",");
        for(int i=1;i<mainModel.yDim; i++) {
            karyotypeDist.Write("\n");
            karyotypeDist.WriteDelimit(outputData[i],",");
            //exportData.exportData(matrix3);
        }
        karyotypeDist.Close();
       */

        //Writing cell locations in (x,y) coordinates and their corresponding karyotype values
        // printMatrix(karyotypeAndLocation);

        karyotype_Location.WriteDelimit(karyotypeAndLocation[0], ",");//Written by column, initialize by writing first column
        for(int i=1;i<(G.xDim*G.yDim); i++) { //loop over remaining columns
            karyotype_Location.Write("\n"); //Required
            karyotype_Location.WriteDelimit(karyotypeAndLocation[i], ",");
        }
        karyotype_Location.Close();


        /**    collecting data on Oxygen conc per grid cell*/
        for(int i=0;i<G.xDim;i++){
            for(int j=0;j<G.yDim;j++){
                O2ConcentrationPerGridCell[i][j]=G.resources.pdegrid2d.Get(i,j);
            }
        }




    }

    public void saveO2Data(int currentTime, OnLatticeGrid G){
        ExportData exportData = new ExportData();
        /** Saves oxygen concentration per grid cell  */
        FileIO OxygenDist=new FileIO(G.mainDir.concat("/output/O2/"+G.tIdx+"_"+currentTime+"_O2Distribution.csv"), "w");


        double [][] O2ConcentrationPerGridCell = new double[G.xDim][G.yDim];
        double [][] TotalO2ConcentrationInTheGrid=new double[(int)G.tEnd][2];//


        /**    collecting data on Oxygen conc per grid cell*/
        for(int i=0;i<G.xDim;i++){
            for(int j=0;j<G.yDim;j++){
                O2ConcentrationPerGridCell[i][j]=G.resources.pdegrid2d.Get(i,j);
            }
        }

        OxygenDist.WriteDelimit(O2ConcentrationPerGridCell[0],",");//
        for(int i=1;i<(G.xDim); i++) { //loop over remaining columns
            OxygenDist.Write("\n"); //Required
            OxygenDist.WriteDelimit( O2ConcentrationPerGridCell[i], ",");
        }
        OxygenDist.Close();
    }







    public void saveMax_deltas(ArrayList<Double[]> max_deltas){
        FileIO DeltasWriter = new FileIO(G.mainDir.concat("/output/deltas.csv"), "w");
        String [] DeltaNames={"Time","max_delta"};



        double[][] doubleArray = new double[max_deltas.size()][2];

        for (int i = 0; i < max_deltas.size(); i++) {
            Double[] temp = max_deltas.get(i);  // Get the array [Time, max_delta]
            doubleArray[i][0] = temp[0];  // Time
            doubleArray[i][1] = temp[1];  // max_delta
        }


        DeltasWriter.WriteDelimit(DeltaNames,",");

        for(int i=0;i< doubleArray.length;i++){
            DeltasWriter.Write("\n");
            String[] row = {Double.toString(doubleArray[i][0]), Double.toString(doubleArray[i][1])};
            DeltasWriter.WriteDelimit(row,",");
        }
        DeltasWriter.Close();

    }
    public void saveMaxDeltasCustom(ArrayList<Double[]> max_deltas){
        FileIO cusDeltasWriter = new FileIO(G.mainDir.concat("/output/max_delta.csv"), "w");
        String [] DeltaNames={"Time","index","x","y","max_delta","Eq_reached"};



        double[][] doubleArray = new double[max_deltas.size()][6];

        for (int i = 0; i < max_deltas.size(); i++) {
            Double[] temp = max_deltas.get(i);  // Get the array [Time, max_delta]
            doubleArray[i][0] = temp[0];  // Time
            doubleArray[i][1] = temp[1];  // max_delta
            doubleArray[i][2] = temp[2];// index
            doubleArray[i][3] = temp[3];
            doubleArray[i][4] = temp[4];
            doubleArray[i][5] = temp[5];

        }


        cusDeltasWriter.WriteDelimit(DeltaNames,",");

        for(int i=0;i< doubleArray.length;i++){
            cusDeltasWriter.Write("\n");
            String[] row = {Double.toString(doubleArray[i][0]),
                    Double.toString(doubleArray[i][1]),
                    Double.toString(doubleArray[i][2]),
                    Double.toString(doubleArray[i][3]),
                    Double.toString(doubleArray[i][4]),
                    Double.toString(doubleArray[i][5])};
            cusDeltasWriter.WriteDelimit(row,",");
        }
        cusDeltasWriter.Close();

    }


    public void savePopulation(double[][] population){
        FileIO Population_output = new FileIO(G.mainDir.concat("/output/PopSize.csv"), "w");
       // FileIO O2Concentration_output = new FileIO(G.mainDir.concat("/output/O2/O2Conc.csv"), "w");
        //FileIO mis_rate_write = new FileIO(G.mainDir.concat("/output/mis_rate.csv"), "w");
        //FileIO div_rate_write = new FileIO(G.mainDir.concat("/output/div_rate.csv"), "w");
        //FileIO OxygenConsumption=new FileIO(G.mainDir.concat("/output/O2/O2Consumption.csv"), "w");


        //int[][] output=new int[population.length][2];
        String [] PopColNames={"Time","N_divRate","C_divRate","N_deathRate","C_deathRate","N_pop","C_pop","Divisions",
                "Deaths","Population","oxygen"};
        Population_output.WriteDelimit(PopColNames,",");

        for(int i=0;i<population.length;i++){
            Population_output.Write("\n");
            Population_output.WriteDelimit(population[i],",");
        }
        Population_output.Close();


/**
        String [] O2ColNames={"Time","O2_conc"};
        O2Concentration_output.WriteDelimit(O2ColNames,",");
        for(int i=0;i<O2Concentration.length;i++){
            O2Concentration_output.Write("\n");
            O2Concentration_output.WriteDelimit(O2Concentration[i],",");

        }
        O2Concentration_output.Close();
*/




        String [] div_rateColNames={"Time","N_divRate","C_divRate","N_deathRate","C_deathRate","N_pop","C_pop"};
        //div_rate_write.WriteDelimit(div_rateColNames,",");



    }


    public void saveO2Consumption(double[][] O2Consumption){

        FileIO OxygenConsumption=new FileIO(G.mainDir.concat("/output/O2/O2Consumption.csv"), "w");

        //int[][] output=new int[population.length][2];
        String [] PopColNames={"Time","Population"};
        OxygenConsumption.WriteDelimit(O2Consumption[0],",");


        for(int i=1;i<O2Consumption.length;i++){
            OxygenConsumption.Write("\n");
            OxygenConsumption.WriteDelimit(O2Consumption[i],",");
        }
        OxygenConsumption.Close();

    }

    //@TODO modify main model to allow image saving inside this method
    public void SaveTumourImage(int currTimeIdx, UIGrid CellGrid, UIGrid pdegrid) {
        OnLatticeGrid G=new OnLatticeGrid(Params.xDim, Params.yDim);
        if (G.imageFrequency > 0 && (currTimeIdx % (int) (G.imageFrequency/G.dt)) == 0) { //Save images at regualar intervals
            CellGrid.ToPNG(this.G.mainDir.concat("/output/img/img_t_"+currTimeIdx*G.dt+".png"));
            System.out.print("Saving image");
        } else if(G.imageFrequency<0 && currTimeIdx ==1 ){ // save image only in the beginning and end
            CellGrid.ToPNG(G.mainDir.concat("/output/img/img_t_"+currTimeIdx*G.dt+".png"));
            pdegrid.ToPNG(G.mainDir.concat("/output/img/diff_img_t_"+currTimeIdx*G.dt+".png"));

        }else if(G.imageFrequency<0 && currTimeIdx ==G.tEnd ){ // save image at the end
            CellGrid.ToPNG(G.mainDir.concat("/output/img/img_t_"+currTimeIdx*G.dt+".png"));
            pdegrid.ToPNG(G.mainDir.concat("/output/img/diff_img_t_"+currTimeIdx*G.dt+".png"));
        }
    }







    public void setWriteDrawValue(int cnv, int location) {

        String temp=(Integer.toString(cnv))+(Integer.toString(location));
        karyotypeDrawValue=Integer.valueOf(temp);
        karyotypeWriteValue=(Integer.toString(cnv))+","+(Integer.toString(location));
    }


    public int[] localRegionIndexes(OnLatticeGrid G, int x,int y, int radius){ //Computes total grid points within a circle centered at (x,y) with radius r
        int Radius=radius;
        int[] localRegion = CircleHood(true, Radius); //indexes for centered circular initialization
        int hoodSize = G.MapHood(localRegion, x, y);
        int[] localRegionIndex = localRegion;
        return localRegionIndex;
    }
    public double[][] distanceDensityData(OnLatticeGrid G){ //Write Distance from vessel, and density of cells in chosen regions
        int[][] center = params.centers;

        int[][] vessels = resources.vessels;

        double[][] DistanceDensity=new double[center.length][6]; //Records: centerX,centerY, average spatial dts from closest vessel,
                                                                 // Density,Diversity at local region, Avg Karyotype dts from diploid
        List<int[][]> CellsAcrossAllSites= new ArrayList<>(); //Holds a list of matrices, each matrix contains cells at a specific site
        List<double[][]> CellToCellData= new ArrayList<>(); //



        double[] AvgKaryoDistanceVector=new double[center.length]; //Stores the average karyotype distance (from diploid) of all cells in EACH chosen region

        for(int k=0;k<center.length;k++) {
            int[] localRegionIndex = localRegionIndexes(G, center[k][0], center[k][1],params.Radius); //Gets all grid cells within the chosen circle of radius R, centered at (x,y)

            // localRegionIndex  contains indexes followed by their (x,y) coordinates.
            //We are only interested in the indexes for now, so we delete x-y coordinates
            int localCellSites = (localRegionIndex.length) / 3;//Contains cell location in the site. The first 1/3 part contains the
                                                               // location indexes and the remaining 2/3 contain (x,y, values)
            int[][] CellsAtSite=new int[localCellSites][24]; //Contains actual cells per site, contains each cell's x,y coordinates followed by karyotype

            int density = 0;
            double[] karyotypeDistances = new double[localCellSites]; //holds the karyotype distances  from diploid of each cell at a given site
            List<Integer> CellPerSiteContainer= new ArrayList<Integer>();
            List<Integer> LocationsWithACell= new ArrayList<Integer>();

            for (int i = 0; i < localCellSites; i++) {

               double agentDistanceFromDiploid = 0;

                if (G.GetAgent(localRegionIndex[i]) != null) {

                    density += 1;
                    Cell agent = G.GetAgent(localRegionIndex[i]);//get the cell at location
                    CellPerSiteContainer.add(agent.hash); //count unique agent hash
                    CellsAtSite[i][0]=G.ItoX(localRegionIndex[i]);
                    CellsAtSite[i][1]=G.ItoY(localRegionIndex[i]);
                    for(int n=2;n<24;n++) {CellsAtSite[i][n]=agent.karyotype[n-2];}

                    LocationsWithACell.add(localRegionIndex[i]); // Record the grid cell
                    agentDistanceFromDiploid = helper.Distance(agent.karyotype,G.diploid_karyotype,"manhattan");//compute agents karyotype distance from diploid


                }

                karyotypeDistances[i] = agentDistanceFromDiploid;


            }
            //helper.printMatrix(helper.deleteZeroRows(CellsAtSite));
            CellsAcrossAllSites.add(helper.deleteZeroRows(CellsAtSite));

            int[] totalCellsAtLocalRegion= helper.toIntArray(CellPerSiteContainer); //Convert integer list to int array. Holds the unique hash codes of cell at each local region
            int DiversityAtRegion= helper.countDistinctElements(totalCellsAtLocalRegion); //Counts the number of distinct cells per the local region

            int[] ActiveCellSites=helper.toIntArray(LocationsWithACell);

            if (density != 0) {
                AvgKaryoDistanceVector[k] = helper.vectorSum(karyotypeDistances) / density; //Find the average karyotype distance of all cells in the region

            } else {
                AvgKaryoDistanceVector[k] = 0;
            }

            double[][] distanceArray = new double[vessels.length][ActiveCellSites.length]; // (size: Vessels X sites) Stores spatial distances of all cells within the current region to a vessel


             for (int i = 0; i < vessels.length; i++){ //loop through all vessels
                for (int j = 0; j < ActiveCellSites.length; j++) { //loop through cells in the region
                    distanceArray[i][j] = G.Dist(G.ItoX(ActiveCellSites[j]), G.ItoY(ActiveCellSites[j]), vessels[i][0], vessels[i][1]); //distance between cell's location and vessel
                }
             }


             double[] AvgDtsOverVessesls=new double[vessels.length];
             if(vessels.length==1){
                AvgDtsOverVessesls[0]=helper.vectorSum(distanceArray[0])/ density;
             }else if(vessels.length>1)
             {for(int i=0;i<(vessels.length);i++){

                 int[][] a={{1,2}};
                 AvgDtsOverVessesls[i]= helper.vectorSum(distanceArray[i])/ density;// Returns mean distance from ith vessel for cells in the region
             }
             }


          double closestVesselDistance = Arrays.stream(AvgDtsOverVessesls).min().getAsDouble(); //Get the minimum of all the distances

            //double distance = Math.min(AvgDtsOverVessesls);
            //double[] densityDistance = {density, distance};

            DistanceDensity[k][0] = center[k][0];
            DistanceDensity[k][1] = center[k][1];
            DistanceDensity[k][2] = closestVesselDistance;
            DistanceDensity[k][3] = density;
            DistanceDensity[k][4] = DiversityAtRegion;
            DistanceDensity[k][5] = AvgKaryoDistanceVector[k];

        }

        int no_of_site_pairs=(center.length-1)*(center.length)/2; //Number of site pairs for n sites is the sum of fist (n-1) whole numbers (Arithmetic sequence)


         double[][] pairwise_Space_karyotye_dts= new double[no_of_site_pairs][3]; // col1 holds  spatial distance between center of two spots,
                                                                                   // col2 holds average spatial distance between cells in spot i and j,
                                                                                    // col 3 and max karyotype distance between cells in spot i and spot j
        //for(int k=0;k<(center.length-1);k++){
        int index=0;
        for(int n=0;n<CellsAcrossAllSites.size();n++) {
            int[][] site_n_locations = Arrays.stream(CellsAcrossAllSites.get(n)) //Extract columns 1 and 2 containing position index of cells at location n
                    .map(row -> Arrays.copyOfRange(row, 0, 2))
                    .toArray(int[][]::new);


            int[][] site_n_karyotypes = helper.getSubmatrix(CellsAcrossAllSites.get(n)); //get corresponding karyotype in columns 3:24


            for (int m = n+1; m < CellsAcrossAllSites.size(); m++) {
                int[][] site_m_karyotypes = helper.getSubmatrix(CellsAcrossAllSites.get(m));

                int[][] site_m_locations = Arrays.stream(CellsAcrossAllSites.get(m)) //Extract columns 1 and 2 containing position index
                        .map(row -> Arrays.copyOfRange(row, 0, 2))
                        .toArray(int[][]::new);



                  double dist_btn_centroids=helper.Distance(center[n],center[m],"euclidean"); // holds the distance between the centers of current spots m and spot n


                double[] Pairwise_spatial_dist = helper.calculatePairwiseDistances(site_n_locations, site_m_locations,"euclidean");
                double[] Pairwise_karyotype_dist = helper.calculatePairwiseDistances(site_n_karyotypes, site_m_karyotypes,"manhattan");
                double[][]pairwise_space_karyo = {Pairwise_spatial_dist,Pairwise_karyotype_dist };
                CellToCellData.add((double[][])helper.Transpose(pairwise_space_karyo)); //Hold pairwise cell-to-cell spatial and karyotype  distances (cell pairs)

                double[]nextRow= {dist_btn_centroids,ArrayMean(Pairwise_spatial_dist),ArrayMax(Pairwise_karyotype_dist)}; //Hold pairwise spatial and karyotype dist (spot pairs)
                pairwise_Space_karyotye_dts[index++] = nextRow;
            }
          }


        //}
        double[][] cell_to_cell_output_data=helper.unlistAndStackMatrices(CellToCellData);
        //System.out.println("cell to cell");
        //helper.printDoubleMatrix(cell_to_cell_output_data);



        String [] colNames={"centerX","centerY", "Avg_Dts_to_vssl","Density","Diversity","Avg_Karyo_dts_from_diploid"};
        FileIO DistDensityWriter=new FileIO(G.mainDir.concat("/"+G.tIdx+"_DistanceDensity.csv"),"w");
        DistDensityWriter.WriteDelimit(colNames , ",");//Written by column, initialize by writing first column
        for(int i=0;i<(DistanceDensity.length); i++) { //loop over remaining columns
            DistDensityWriter.Write("\n"); //Required
            DistDensityWriter.WriteDelimit(DistanceDensity[i], ",");
        }

        DistDensityWriter.Close();

        String [] columnNames={"Dist_between_centers","Avg_Dts_btn_sites","Max_karyotype_Dts"};
        FileIO SiteSpatialKaryotyppeWriter=new FileIO(G.mainDir.concat("/"+G.tIdx+"_Spatial_karyo_dts.csv"),"w");
        SiteSpatialKaryotyppeWriter.WriteDelimit(columnNames , ",");//Written by column, initialize by writing first column
        for(int i=0;i<(pairwise_Space_karyotye_dts.length); i++) { //loop over remaining columns
            SiteSpatialKaryotyppeWriter.Write("\n"); //Required
            SiteSpatialKaryotyppeWriter.WriteDelimit(pairwise_Space_karyotye_dts[i], ",");
        }
        SiteSpatialKaryotyppeWriter.Close();

        //Record pairwise cell-to-cell spatial and  karyotype distances for all cells in each pair of spots
        String[] cols={"Dist_btn_cells","karyotype_Dts"};
        FileIO cell2cellWriter=new FileIO(G.mainDir.concat("/"+G.tIdx+"_cell2cell.csv"),"w");
        cell2cellWriter.WriteDelimit(cols , ",");//Written by column, initialize by writing first column
        for(int i=0;i<(cell_to_cell_output_data.length); i++) { //loop over remaining columns
            cell2cellWriter.Write("\n"); //Required
            cell2cellWriter.WriteDelimit(cell_to_cell_output_data[i], ",");
        }
        cell2cellWriter.Close();


        return DistanceDensity;

    }



    public void DistanceDensity2(OnLatticeGrid G) {
       // This method writes out output data on distance from vessel, density, diversity etc to be analyzed
        int[] center= {100,100};
        int[]Radii={1};
        int Radius=1;
        int[] localRegion =localRegionIndexes(G,center[0],center[1],Radius) ; //indexes for centered circular initialization
        double[][] DistanceDensity=new double[center.length][6]; //Records: centerX,centerY, average spatial dts from closest vessel,
        // Density,Diversity at local region, Avg Karyotype dts from diploid


        List<int[][]> CellsAcrossAllSites= new ArrayList<>(); //Holds a list of matrices, each matrix contains cells at a specific site
        List<double[][]> CellToCellData= new ArrayList<>(); //

    }

    public void WriteRnNumbers(ArrayList arrayList,ArrayList arrayList2) {

        double[] doubleArray = new double[arrayList.size()];
        double[] doubleArray2 = new double[arrayList2.size()];

        for(int i = 0; i < arrayList.size(); i++) {
            doubleArray[i] = (double) arrayList.get(i);
        }

        for(int i = 0; i < arrayList2.size(); i++) {
            doubleArray2[i] = (double) arrayList2.get(i);
        }

        FileIO RnWriter=new FileIO(G.mainDir.concat("/output/random.csv"),"w");
        FileIO RnWriter2=new FileIO(G.mainDir.concat("/output/Mrandom.csv"),"w");
/**
        RnWriter.WriteDelimit(Collections.singletonList(doubleArray[0]),",");
        for(int i=1;i<doubleArray.length;i++){
            RnWriter.Write("\n");
            RnWriter.WriteDelimit(Collections.singletonList(doubleArray[i]),",");
        }

        RnWriter.Close();
*/

        RnWriter2.WriteDelimit(Collections.singletonList(doubleArray2[0]),",");
        for(int i=1;i<doubleArray2.length;i++){
            RnWriter2.Write("\n");
            RnWriter2.WriteDelimit(Collections.singletonList(doubleArray2[i]),",");
        }

        RnWriter2.Close();

        }









    public void cellTypePerRegionWriter(OnLatticeGrid G) {
        int[][] centers = params.centers;
        for (int k = 0; k < centers.length; k++) {
            cellTypePerRegion(G,centers[k][0],centers[k][1]);
        }

    }

    public void cellTypePerRegion(OnLatticeGrid G, int centerX, int centerY) { //Collect cell data within select regions, records cell specific location and karyotype

      int[] localRegionIndex = localRegionIndexes(G, centerX,centerY,params.Radius);

        // localRegionIndex  containing indexes followed by their (x,y) coordinates.
        //We are only interested in the indexes for now, so we delete x-y coordinates
        int n = (localRegionIndex.length) / 3;//The first 1/3 part contains the location indexes and the remaining 2/3 contain (x,y, values)
        double[][] cellsPerRegion = new double[n][25]; // contains the ff information:  x-coord,y-corrd,distance from vessel, karyotype
        for(int i=0;i<n;i++){
            cellsPerRegion[i][0]=G.ItoX(localRegionIndex[i]); //set x coordinate
            cellsPerRegion[i][1]=G.ItoY(localRegionIndex[i]); //set y cordinate
            Cell agent=G.GetAgent(localRegionIndex[i]); //get cell
                //@TODO modify to work for multyple blood vesses
            cellsPerRegion[i][2]=G.Dist(G.ItoX(localRegionIndex[i]),G.ItoY(localRegionIndex[i]),101, 101);//record spatial distance between agents location and vessel


            if (agent != null) { //if there is an agent, record karyotype
                for(int j=3;j<25;j++){
                    cellsPerRegion[i][j]=agent.karyotype[j-3]; //set the corresponding karyotype
                }
            }else{ //IF there is no agent record zero copy numbers
                for(int j=2;j<24;j++){
                    cellsPerRegion[i][j]=0;
                }

            }

        }
        FileIO CellTypeWriter=new FileIO(G.mainDir.concat("/output/cellType/"+centerX+"_"+centerY+"-"+"t_"+G.tIdx+"_location.csv"),"w");
        String [] columNames={"x","y","vessel dts","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22"};
        CellTypeWriter.WriteDelimit(columNames, ",");//Written by column, initialize by writing first column
        for(int i=0;i<(cellsPerRegion.length); i++) { //loop over remaining columns
            CellTypeWriter.Write("\n"); //Required
            CellTypeWriter.WriteDelimit(cellsPerRegion[i], ",");
        }

        CellTypeWriter.Close();

 }





}