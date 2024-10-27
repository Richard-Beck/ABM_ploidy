package onLatticeCA_jar;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.*;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;


import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.Random;
//import MelanomaModel.Pars;
//import MelanomaModel.OnLatticeVis;

import java.io.File;
import java.util.Arrays;

import java.io.*;

import static HAL.Util.*;

public class OnLatticeGrid extends AgentGrid2D<Cell> implements SerializableModel {


    Params Params = new Params();


    int TotalRun=0;
    ArrayList<Double[]> max_deltas = new ArrayList<>();


    boolean homeostasisReached = false;


    /* ========================================================================
     * --- Parameters ----
     * ======================================================================== */
    public int xDim = Params.xDim; // x-dimension of domain (grid size )
    public int yDim = Params.yDim; // y-dimension of domain (grid size)

    Resources resources;
    int[] initPopSizeArr;
    double initialSize = Params.initialTumorSize; // Initial population size
   static int InitialPopulationType;

    double initialSizeProp; // Initial cell density relative to (physical) carrying capacity
    double rFrac = 0.05; // Initial resistance fraction in [0,1]
    double divisionRate_S = Params.natural_growth_rate; // Proliferation rate of sensitive cells in d^-1. Proliferation will be attempted at this rate.
    double divisionRate_R = Params.natural_growth_rate; // Proliferation rate of resistant cells in d^-1
    double movementRate_S = Params.movementRate_S; // Movement rate of sensitive cells in d^-1. Movement will be attempted at this rate.
    double movementRate_R = 0; // Movement rate of resistant cells in d^-1.
    double deathRate_S = Params.natural_death_rate; // Natural death rate of sensitive cells in d^-1. Mean life span of a cell will be 1/k
    double deathRate_R = Params.natural_death_rate;// Natural death rate of resistant cells in d^-1.
    double drugKillProportion = Params.drugKillProportion;//0.75; // Drug induced death rate in d^-1.
    int[] hood = Util.VonNeumannHood(false); // Define division in von Neumann neighbourhood.
    // Argument indicates that we don't want the center position returned (can only divide into neighbouring squares).
    int tEnd = Params.tEnd;
    /**
     * Originally 7 End time in d
     */

    double dt = Params.dt; // Time step in d
    int nTSteps; // Number of time steps to simulate
    double[][] treatmentScheduleList = {{0, tEnd, 0}}; // Treatment schedule in format {{tStart, tEnd, drugConcentration}}
    int nReplicates = 1; // Number of replicates


    //double[]  consumption=new double[3];

    static String mainDir = "/Users/4477116/Documents/projects/ABM_ploidy"; //Main project directory
    ArrayList<String> randomNumbers = new ArrayList<String>();


    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    int RegularDivisionCounter;
    int missegregationCounter;
    int normalCellCounter,cancerCellCounter;
    int DeadCellCounter;
    int divisionCounter = 0;

    public int monoploidCounter;

    public int diploidCounter;

    helperMethods helper = new helperMethods();


    int[][] baseKaryotype = new int[Params.xDim * Params.yDim][22];
    int[] diploidKaryotype = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};



    /**
     * Select location to initialize cell population [0-6]:
     * 0 = User input; 1= Center, 2=Top right, 3=Top left, 4=Bottom right etc
     */

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // Helper variables
    int tIdx = 0;
    double diffIdx=0.0;
    double currDrugConcentration; // Current drug concentration in [0,1].
    int[] cellCountsArr = new int[2]; // Array to hold number of each cell type

    Rand rn; // Random number generator
    Rand rn_ICs; // Random number generator for initial conditions
    int seed = 42;//Math.abs(random.nextInt());

    int verboseLevel = 1;
    double printFrequency = 3; // Frequency at which output is printed to the screen

    FileIO cellCountLogFile = null;
    double logCellCountFrequency = -1; // Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging
    double[] extraSimulationInfo;
    String[] extraSimulationInfoNames;

    // Output - Visualisation
    UIGrid vis;

    Boolean visualiseB = true; // Whether or not to show visualization
    Boolean visualiseA = true;
    int scaleFactor = 2; // Scale factor used to display grid
    final static int BLACK = Util.RGB(0, 0, 0);
    int imageFrequency = 1; // Frequency at which an image of the tumour is saved. Negative number turns it off

    int[] diploid_karyotype = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    ExportData exporter;
    public ArrayList<Double> mainRnNumbers = new ArrayList<>();

    public static void main(String[] args) throws IOException {
        if(args.length<2){
            System.out.println("intended usage: OnLatticeGrid path2Params outputPath");
        }
        String path2Params = args[0];
        String outputPath = args[1];

        Params Params = new Params();
        Params.load(path2Params);


        OnLatticeGrid myModel = new OnLatticeGrid(Params.xDim, Params.yDim);
        myModel.SetSeed(myModel.seed);
        myModel.exporter = new ExportData(outputPath,path2Params);
        myModel.resources = new Resources(myModel, myModel.xDim, myModel.yDim);;
        myModel.resources.initialize_vessels(Params.path2Vessels);
        myModel.resources.vessels2Indices();

        myModel.Run();
        myModel.Close();
    }


    public OnLatticeGrid(int xDim, int yDim) {
        super(xDim, yDim, Cell.class);
    }




    // Function used as part of SerializableModel to allow saving the model's state so that I can restart
    // the simulation exactly where I left off at the end of the last treatment interval.
    @Override
    public void SetupConstructors() {
        _PassAgentConstructor(Cell.class);
    }



    public void SetSeed(int seed) {
        this.rn = new Rand(seed);
        this.rn_ICs = new Rand(seed);

    }


    public double[] GetModelState() {
        return new double[]{tIdx, tIdx * dt, cellCountsArr[0], cellCountsArr[1], Util.ArraySum(cellCountsArr),
                currDrugConcentration, divisionRate_S, divisionRate_R, movementRate_S, movementRate_R, deathRate_S, deathRate_R, drugKillProportion, dt};

    }

    // ------------------------------------------------------------------------------------------------------------
    // Seeding functions for the cells
    public void InitCells(int N){
        double seedingDensity = (double) (N) /(xDim*yDim);
        for(int i = 0; i<xDim*yDim; i++){
            if(rn.Double()<seedingDensity){
                Cell c = NewAgentSQ(i);
                c.init(0,diploid_karyotype);
                c.Draw();
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    public void StepCells() {
        normalCellCounter=0; cancerCellCounter=0;
        ShuffleAgents(rn); //shuffle order of for loop iteration over cells
        for (Cell c : this) { //iterate over all cells in the grid
            if(c.cancerous) {
                cancerCellCounter++;
            } else normalCellCounter++;
            //Check if cell's chromosamal copy numbers are  within the right range
            if (helper.isAboveThreshold(5, c.karyotype)) { //if a copy number  is above 5, kill cell
                vis.SetPix(c.Isq(), BLACK); // Death
                c.Dispose();
                System.out.println("cell is dead");
                DeadCellCounter += 1;

            } else  c.DivideOrDie();
        }

    }

    // ------------------------------------------------------------------------------------------------------------
    public void Run() {

        // Initialise visualisation window
        this.vis = new UIGrid(xDim, yDim, scaleFactor, visualiseB);
        visualize visuals = null;
        if(Params.runGUI ==1 ) visuals =new visualize(this,"visuals");


        InitCells((int)Params.initialTumorSize);
        homeostasisTests h = new homeostasisTests();
        // Run the simulation
        boolean completedSimulation = false;
        while (!completedSimulation) {
            PrintStatus(tIdx);
            // NOTE:
            //  resources.maximum_delta will be NaN if the O2 initial condition is 0.
            //  This has the potential to cause problems /RJB
            double diffusionTime = (double) tIdx;
            int cxx = 0;
            double AvgO2 = 0.;

            while (Params.diffusion_tol < resources.maximum_delta | cxx <10) {
                cxx++;
                diffusionTime+=resources.setDirechletCond();
                AvgO2=  resources.pdegrid2d.GetAvg();
                exporter.writeOxygenSummary(diffusionTime,AvgO2,resources.maximum_delta);

            }
            resources.Drawvessels();

            if(tIdx == 1){

                exporter.saveImage(vis,resources.currV,tIdx);
            } else if (Params.imageFrequency > 0 && (tIdx % (int) (Params.imageFrequency/dt)) == 0){
                exporter.saveImage(vis,resources.currV,tIdx);
            }

            if (tIdx == 1) {
                if (Params.writeKaryotypeData) {
                    exporter.saveKaryotypeLocationData(tIdx, this);
                    exporter.saveOxygenField(tIdx,resources.pdegrid2d);
                }
            } else if (verboseLevel > 0 && (tIdx % (int) (Params.writeFrequency / dt)) == 0) {
                if (Params.writeKaryotypeData) {
                    exporter.saveKaryotypeLocationData(tIdx, this);
                    exporter.saveOxygenField(tIdx,resources.pdegrid2d);
                }
            }
            if(tIdx%5 == 0&tIdx>1) {
                exporter.createGifFiles(vis,resources.currV);
            }


            if(!homeostasisReached){
                homeostasisReached = h.test(AvgO2, normalCellCounter);
                if(homeostasisReached){ System.out.println("Homeostasis reached");
                }
            }

            tIdx++;
            // Check if the stopping condition is met
            completedSimulation = (tIdx > Params.tEnd);//?true:false;
            StepCells();
            exporter.writeSummary((double) tIdx, normalCellCounter, cancerCellCounter);
        }


        exporter.close();

        if(visuals != null) visuals.close();
    }


    public Boolean SaveCurrentCellCount(int currTimeIdx) {
        Boolean successfulLog = false;
            if ( logCellCountFrequency > 0 &&(currTimeIdx % (int) (logCellCountFrequency / dt)) == 0) {
                cellCountLogFile.WriteDelimit(GetModelState(), ",");
                if (extraSimulationInfoNames != null) {
                    cellCountLogFile.Write(",");
                    cellCountLogFile.WriteDelimit(extraSimulationInfo, ",");
                }
                cellCountLogFile.Write("\n");
                successfulLog = true;
            }
        return successfulLog;
    }

    public void PrintStatus(int currTimeIdx) {

        if (verboseLevel > 0 && (currTimeIdx %  (int) (printFrequency/dt)) == 0) {
            System.out.println("Time: " + currTimeIdx*dt+ " - Population Size: " + Util.ArraySum(cellCountsArr)+
                  " -  Dead Cells: " +  DeadCellCounter + " - missegregations: " + missegregationCounter+
                          " - Regular Divisions: " + RegularDivisionCounter+ " - Mono: "+ monoploidCounter+ " - Diploid: " + diploidCounter
            +"  -Total population: "+Pop() );


        }
    }

    // ----------------------------------------------------------------------------------------------------------
    // Function to clean up loose ends after running a simulation.
    public void Close() {
        if (cellCountLogFile!=null) cellCountLogFile.Close();

    }

    public void Close(Boolean logged) {
        if (!logged) {
            tIdx--;
            SaveCurrentCellCount(0);}
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }


}