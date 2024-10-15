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
    Params par = new Params();

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

    static String mainDir = "C:/Users/4473331/Documents/projects/008_birthrateLandscape/ABM_ploidy"; //Main project directory
    ArrayList<String> randomNumbers = new ArrayList<String>();


    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    int RegularDivisionCounter;
    int missegregationCounter;

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
    Random random = new Random();
    int seed = 42;//Math.abs(random.nextInt());

    int verboseLevel = 1;
    double printFrequency = 3; // Frequency at which output is printed to the screen
    double writeFrequency = 5;//Params.tEnd; //Frequency at which output is written to csv file
    String cellCountLogFileName;
    FileIO cellCountLogFile = null;
    double[][] outputArr; // Array to hold the simulation results if they're not directly logged to file
    double logCellCountFrequency = -1; // Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging
    double[] extraSimulationInfo;
    String[] extraSimulationInfoNames;

    // Output - Visualisation
    UIGrid vis;
    UIGrid vis2;


    Boolean visualiseB = true; // Whether or not to show visualization
    Boolean visualiseA = true;
    int scaleFactor = 2; // Scale factor used to display grid
    int pause = 2; // 50  Pause between time steps to smoothen simulation
    final static int BLACK = Util.RGB(0, 0, 0);
    final static int RED = Util.RGB(1, 0, 0);
    final static int WHITE = Util.RGB(1, 1, 1);
    int imageFrequency = 1; // Frequency at which an image of the tumour is saved. Negative number turns it off

    String imageOutDir = "C:/Users/4473331/Documents/projects/008_birthrateLandscape/ABM_ploidy"; // Directory which to save images to

    // Output - Model file
    Boolean saveModelState = false;
    Boolean fromScratch = false;
    String savedModelFileName = null; // Name of model file to load when continuing a previous run
    int[] diploid_karyotype = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    int[][] inputKaryotype = Params.inputKaryotype;

    ExportData Export = new ExportData();
    ExportData exporter;
    public ArrayList<Double> CellRnNumbers = new ArrayList<>();
    public ArrayList<Double> mainRnNumbers = new ArrayList<>();


    // ------------------------------------------------------------------------------------------------------------



    public static void main(String[] args) throws IOException {
        if(args.length<2){
            System.out.println("intended usage: OnLatticeGrid path2Params outputPath");
        }
        String path2Params = args[0]; //"C:/Users/4473331/Documents/projects/008_birthrateLandscape/ABM_ploidy/input/parameters.txt";
        String outputPath = args[1]; //"C:/Users/4473331/Documents/projects/008_birthrateLandscape/ABM_ploidy/output";

        Params Params = new Params();
        Params.load(path2Params);
        //Params.readParameters(path2Params);



        InitialPopulationType=Params.InitialPopulationType;

        Params.readVessel();



        /** Begin loop here when running replicates with various parameter values. Example: replicates with diff misSegRates */
        // for(int rat=0;rat<Params.GrowthDeathRates.length;rat++){
       // for( double conrat:Params.consumptionrates){
        //    Params.oxygen_consumption_rate=conrat;
        //     Params.natural_growth_rate=Params.GrowthDeathRates[rat][0];
        //     Params. natural_death_rate=Params.GrowthDeathRates[rat][1];
        //for(int pop=0;pop<Params.InitialSizes.length;pop++){
        //  Params.initialTumorSize=Params.InitialSizes[pop];
        //for (double tend : Params.tEndValues) {
       // for (int  vesselOxygen: Params.oxygenAtVessels) {
        //for(int oxygen_level:Params.initial_levels){
        //  Params.initial_oxygen_conc=oxygen_level;
         // Params.vssl_bdary_value = vesselOxygen;
        /** Begin loop here when running replicates with various parameter values*/


        OnLatticeGrid myModel = new OnLatticeGrid(Params.xDim, Params.yDim);
        myModel.exporter = new ExportData(outputPath,path2Params);


        double[][] mis_rate = new double[myModel.tEnd][2];
        helperMethods.readAndAppend(myModel.seed, myModel.mainDir.concat("/input/seed.txt"));

        // Runs an example simulation using the default parameters
        // Setup


        int nReplicates = 1;
        int[] replicateIdList = new int[]{myModel.seed};
        int replicateId;
        int[] initPopSizeArr = {0, 0};
        double[][] treatmentScheduleList = {{0, myModel.tEnd, 0}}; // Treatment schedule in format {{tStart, tEnd, drugConcentration}}
        String outDir = "./tmp/";
        Boolean fromScratch = true;
        Boolean visualiseB = true;
        Boolean visualiseA = true;
        //int pause = 100; //250;
        String savedModelFileName = null; // Name of model file to load when continuing a previous run
        String currSavedModelFileName;

        // Update model object with new params
        myModel.fromScratch = fromScratch;


        // Run the sweep
        System.out.println("Number of replicates " + nReplicates);
        if (nReplicates > 1) {
            replicateIdList = new int[nReplicates];
            for (int i = 0; i < nReplicates; i++) {
                replicateIdList[i] = i;
            }
        }
        if (myModel.initialSizeProp > 0) { // Make it so that can use proportional definition of initial density
            myModel.initialSize = (int) (myModel.initialSizeProp * myModel.xDim * myModel.yDim);
        }
        initPopSizeArr = new int[]{(int) Math.floor(myModel.initialSize * (1 - myModel.rFrac)), (int) Math.ceil(myModel.initialSize * myModel.rFrac)};
        //new File(outDir).mkdirs();
        for (int replicateIdx = 0; replicateIdx < nReplicates; replicateIdx++) {
            replicateId = replicateIdList[replicateIdx];
            currSavedModelFileName = savedModelFileName == null ? outDir + "RepId_" + replicateId + ".bin" : savedModelFileName;


            // Models can be created either from scratch, or be loaded from a previous run
            if (fromScratch) {
                // Set the random number seed
                if (nReplicates == 1) {

                    myModel.SetSeed(myModel.seed);
                } else {

                    myModel.SetSeed(replicateId, myModel.seed);
                }

                // Set the logging behaviour
                // myModel.ConfigureVisualisation(visualiseB, pause);
                if (myModel.imageOutDir.length() > 0) {
                    myModel.visualiseB = true;
                    myModel.visualiseA = true;
                    String currImgOutDir = outDir + "/img/" + "RepId_" + replicateId;
                    //new File(currImgOutDir).mkdirs();
                    myModel.ConfigureImaging(currImgOutDir, myModel.imageFrequency);
                }


                // Initialise the simulation
               // myModel.InitialiseCellLog(outDir + "RepId_" + replicateId + ".csv");
                myModel.SetInitialState(initPopSizeArr);
            } else {
                myModel = LoadState(currSavedModelFileName);
            }

            // Run the simulation
            myModel.SetTreatmentSchedule(treatmentScheduleList);
            myModel.Run();
            myModel.Close();
            if (myModel.saveModelState) {
                myModel.SaveModelState(currSavedModelFileName);
            }
        }


    /**  End loop here when running replicated for various parameter values*/
   // }

    /**
     * End loop here when running replicated for various parameter values
     */

    }
    // ------------------------------------------------------------------------------------------------------------
    // Grid Constructors



    public OnLatticeGrid(int xDim, int yDim) {
        super(xDim, yDim, Cell.class);
        resources = new Resources(this, xDim, yDim);


    }

    public OnLatticeGrid(int x, int y, double[] paramArr, double dt) {
        super(x, y, Cell.class);
        SetParameters(paramArr);
        this.dt = Params.dt;


    }


    // Function used as part of SerializableModel to allow saving the model's state so that I can restart
    // the simulation exactly where I left off at the end of the last treatment interval.
    @Override
    public void SetupConstructors() {
        _PassAgentConstructor(Cell.class);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Functions to parameterise the model
    public void SetParameters(double[] paramArr) {
        this.divisionRate_S = paramArr[0];
        this.divisionRate_R = paramArr[1];
        this.movementRate_S = paramArr[2];
        this.movementRate_R = paramArr[3];
        this.deathRate_S = paramArr[4];
        this.deathRate_R = paramArr[5];
        this.drugKillProportion = paramArr[6];
        //this.karyotype=paramArr[7];

    }

    public void SetInitialState(int[] initialStateArr) {
        this.initPopSizeArr = initialStateArr;
    }

    public void SetDrugConcentration(double drugConcentration) {
        this.currDrugConcentration = drugConcentration;
    }

    public void SetTreatmentSchedule(double[][] treatmentScheduleList) {
        this.treatmentScheduleList = treatmentScheduleList;
    }

    public void SetSeed(int seed) {
        this.rn = new Rand(seed);
        this.rn_ICs = new Rand(seed);

    }

    public void SetSeed(int seed_Simulation, int seed_ICs) {
        this.rn = new Rand(seed_Simulation);
        this.rn_ICs = new Rand(seed_ICs);
    }

    public void SetVerboseness(int verboseLevel, double printFrequency){
        this.verboseLevel = verboseLevel;
        this.printFrequency = printFrequency;
    }

    public void ConfigureVisualisation(boolean visualiseB, int pause) {
        this.visualiseB = visualiseB;
        this.visualiseA = visualiseA;
        this.pause = pause;
    }

    public void ConfigureImaging(String imageOutDir, int imageFrequency) {
        /*
         * Configure location of and frequency at which tumour is imaged.
         */
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
    }

    public double[] GetParameters() {
        return new double[]{divisionRate_S, divisionRate_R, movementRate_S, movementRate_R,
                deathRate_S, deathRate_R, drugKillProportion};
    }

    public double[] GetModelState() {
        return new double[]{tIdx, tIdx * dt, cellCountsArr[0], cellCountsArr[1], Util.ArraySum(cellCountsArr),
                currDrugConcentration, divisionRate_S, divisionRate_R, movementRate_S, movementRate_R, deathRate_S, deathRate_R, drugKillProportion, dt};

    }

    // ------------------------------------------------------------------------------------------------------------
    // Seeding functions for the cells

    public void InitSimulation_Specific(int s0, int r0) {
        OnLatticeGrid onlatticeGrid = new OnLatticeGrid(xDim, yDim);
        Cell cell = new Cell();

        Initialization initialization = new Initialization(this);

        try {

            resources.InitialO2Concentrations(); //Set initial oxygen concentration first. Required to compute division Rates on initial population
        } catch (Exception e) {
            throw new RuntimeException(e);
        }


        //Fill basekaryootype with 2's
        for (int i = 0; i < baseKaryotype.length; i++) {
            baseKaryotype[i] = diploidKaryotype;
        }


        System.out.println("Initial simulation");
        //places tumor cells randomly on the dish
        int[] distributionOfResistanceArr = new int[s0 + r0];
        int[] initialLocArr = new int[s0 + r0];// int[xDim * yDim];
        //int[]  initialLocArro=new int[s0+r0];
        Arrays.fill(cellCountsArr, 0); // clear the cell counter


        // Generate a list of random assignment to sensitive or resistance
        for (int i = 0; i < s0; i++) {
            distributionOfResistanceArr[i] = 0;
        }
        for (int i = s0; i < s0 + r0; i++) {
            distributionOfResistanceArr[i] = 1;
        }
        int[] distArr = distributionOfResistanceArr;


        rn_ICs.Shuffle(distributionOfResistanceArr);


       // helperMethods.shuffleArray(distributionOfResistanceArr, rn);




        int[][] InitialKaryotype = baseKaryotype;
        if (InitialPopulationType == 0) {

            InitialKaryotype = initialization.getIniPopKaryotype();

        } else {
            InitialKaryotype = baseKaryotype;

        }


        int[] initialPopLocArr = initialization.setLocation();
        System.out.println(" locations read " + initialPopLocArr.length);

        double Intended_radius = Math.ceil((double) Math.sqrt(Params.initialTumorSize / Math.PI));
        int gridRadius = Math.min(Params.xDim / 2, Params.yDim / 2);
        double correction = Math.ceil(Math.max((Intended_radius - gridRadius), 0));
        //Generate initial cell population
        for (int i = 0; i < s0 + r0; i++) {
            // replace tumorNeighborhood with initialLocArro to initialize with population in a square
            Cell c = NewAgentSQ(initialPopLocArr[i]); // NewAgentSQ(viableInitialLocArro[i]);
            c.init(distributionOfResistanceArr[i],
                    InitialKaryotype[i]);

            c.Draw();
            cellCountsArr[c.resistance] += 1; // Update population counter


        }
        if (Params.writeCellTypeData) {
            Export.cellTypePerRegionWriter(this);
        }
         //Record cell and vessel locations in select subregions (of original data)
        // for comparison with abm output later
        // Export.distanceDensityData(this); //Write data on vessel dts and cell density in select regions.

        diploidCounter = (r0 + s0);


    }


    // ------------------------------------------------------------------------------------------------------------
    public void StepCells() {
        Arrays.fill(cellCountsArr, 0);//clear the cell counts
        ShuffleAgents(rn); //shuffle order of for loop iteration over cells
        for (Cell c : this) { //iterate over all cells in the grid
            cellCountsArr[c.resistance] += 1; // Update population counter
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

        try {
            resources.InitVessels();
        } catch (IOException e) {
            throw new RuntimeException(e);

        }


       //GifMaker myGif1=new GifMaker(mainDir.concat("/output/Gif/"+Params.initialTumorSize+"_Cells.gif"),1,true);
        //GifMaker myGif2=new GifMaker(mainDir.concat("/output/Gif/"+Params.initialTumorSize+"_Vessels.gif"),1,true);
        System.out.println("Time "+tIdx+"__This is dt "+dt+" --and size "+Math.floor((tEnd*nReplicates)/dt));
        double[][] Population=new double[(int) (Math.floor((tEnd*nReplicates)/dt)+1)][11]; //Record population at each time point


       // ExportData Export = new ExportData();



        int cancerPop=0;
        int normalPop=(int)initialSize;




        // Initialise visualisation window
         this.vis = new UIGrid(xDim, yDim, scaleFactor, visualiseB);


        visualize visuals = null;
        if(Params.runGUI ==1 ) visuals =new visualize(this,"visuals");


        Boolean completedSimulationB = false;
        Boolean logged = false;
        currDrugConcentration = treatmentScheduleList[0][2];
        double totalO2;
        double AvgO2;

        // Set up the grid and initialise log if this is the beginning of the simulation
        if (tIdx == 0) {

            System.out.println("Begining delta "+resources.pdegrid2d.MaxDelta());

            // Set vessel locations
            System.out.println("number of vessels: "+resources.vessels.length);

            Population[0][0]=tIdx;
            Population[0][1]=Params.natural_growth_rate;
            Population[0][2]=Params.natural_growth_rate;
            Population[0][3]=Params.natural_death_rate;
            Population[0][4]=Params.natural_death_rate;
            Population[0][5]=initialSize;
            Population[0][6]=0;
            Population[0][7]=0;
            Population[0][8]=0;
            Population[0][9]=(int)initialSize;
            Population[0][10]=xDim*yDim*Params.initial_oxygen_conc;




            InitSimulation_Specific(initPopSizeArr[0], initPopSizeArr[1]);
            PrintStatus(0);


           // if (cellCountLogFile == null && cellCountLogFileName != null) {
             //   InitialiseCellLog(this.cellCountLogFileName);
            //}
           // SaveCurrentCellCount(0);
            //SaveTumourImage(tIdx);
            tIdx = 1;

        } else {
            // Continue from a restart
            if (cellCountLogFileName != null) {
                cellCountLogFile = new FileIO(cellCountLogFileName, "a");
            }
            for (Cell c : this) {
                c.Draw();
            }
        }


        // Run the simulation
        double currIntervalEnd;
        int initialCellNumber = Util.ArraySum(cellCountsArr);
        if (treatmentScheduleList == null) treatmentScheduleList = new double[][]{{0, tEnd, currDrugConcentration}};

        for (int intervalIdx = 0; intervalIdx < treatmentScheduleList.length; intervalIdx++) {
            currIntervalEnd = treatmentScheduleList[intervalIdx][1];
            nTSteps = (int) Math.ceil(currIntervalEnd / dt);

            currDrugConcentration = treatmentScheduleList[intervalIdx][2];
            completedSimulationB = false;
            homeostasisTests h = new homeostasisTests();

            while (!completedSimulationB) {
                totalO2 = helper.getTotalO2ConcInTheGrid(resources.pdegrid2d);
                AvgO2 = helper.getTotalO2ConcInTheGrid(resources.pdegrid2d) / 40000;
                //if(tIdx==9){System.out.println("RN3 "+rn.Double());}

                vis.TickPause(pause);

                PrintStatus(tIdx);




                double cancerCellDivRate = 0;
                double normalCellDivRate = 0;
                double normalCellmis_Rate = 0;
                double cancerCellmis_Rate = 0;
                int CancerCellCounter = 0;
                int NormalCellCounter = 0;
                double NormalCellDeath = 0;
                double CancerCellDeath = 0;


                if ((tIdx > Params.Resource_stop)) {
                    Params.Resource_on = false;
                } else {
                    Params.Resource_on = true;
                }



                logged = SaveCurrentCellCount(tIdx);






                // NOTE:
                //  resources.maximum_delta will be NaN if the O2 initial condition is 0.
                //  This has the potential to cause problems /RJB
                double diffusionTime = (double) tIdx;
                int cxx = 0;
                while (Params.diffusion_tol < resources.maximum_delta | cxx <10 | !Double.isFinite(resources.maximum_delta)) {
                    cxx++;
                    diffusionTime+=resources.setDirechletCond();
                    AvgO2=helper.getTotalO2ConcInTheGrid(resources.pdegrid2d)/(xDim*yDim);
                    exporter.writeOxygenSummary(diffusionTime,AvgO2,resources.maximum_delta);
                    //System.out.println(resources.maximum_delta);
                    diffIdx = helper.roundToDecimalPlaces((diffIdx +Params.diffusion_dt),2);
                }
                resources.Drawvessels();
                for (Cell agent : this) {
                    if (agent.cancerous) {
                        cancerCellmis_Rate += agent.missegregationRate;
                        cancerCellDivRate += agent.divisionRate;
                        CancerCellCounter += 1;
                        CancerCellDeath += agent.deathRate;

                    } else {
                        normalCellDivRate += agent.divisionRate;
                        normalCellmis_Rate += agent.missegregationRate;
                        NormalCellCounter += 1;
                        NormalCellDeath += agent.deathRate;
                    }
                }

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


                TotalRun = 0;

                if (TotalRun > 200 && resources.maximum_delta > 5) {
                    System.out.println("Equilibrium not reached"); // proceed regardless
                }


                if(verboseLevel > 0 && (tIdx % 1)== 0){
                    Population[tIdx][0]=tIdx;
                    Population[tIdx][1]=normalCellDivRate/ NormalCellCounter; //Effective div. rate(sum of all normal cell div. rates divided total Normal cells )
                    Population[tIdx][2]=cancerCellDivRate/ CancerCellCounter;//Effective div. rate(sum of all cancer cell div. rates divided total cancer cells)
                    Population[tIdx][3]=NormalCellDeath/NormalCellCounter; //Normal cell Effective death rate
                    Population[tIdx][4]=CancerCellDeath/CancerCellCounter;//Cancer cell Effective death rate
                    Population[tIdx][5]=NormalCellCounter;       //Total normal cells
                    Population[tIdx][6]=CancerCellCounter; //Total cancer cells
                    Population[tIdx][7]=divisionCounter; //Total Divisions
                    Population[tIdx][8]=DeadCellCounter; //Total dead
                    Population[tIdx][9]=this.Pop();     //Total population
                    Population[tIdx][10]=helper.getTotalO2ConcInTheGrid(resources.pdegrid2d); //Total oxygen
                }


                 cancerPop +=CancerCellCounter;
                 normalPop += NormalCellCounter;


                homeostasisReached = h.test(AvgO2,NormalCellCounter);

                tIdx++;
                // Check if the stopping condition is met
                completedSimulationB = (tIdx > nTSteps);//?true:false;
                StepCells();
                exporter.writeSummary((double) tIdx, NormalCellCounter, CancerCellCounter);
            }

            //myGif1.Close();
            //myGif2.Close();
            //Export.savePopulation(Population);
            //Export.saveMax_deltas(max_deltas);

           // if(Params.writeCellTypeData){Export.cellTypePerRegionWriter(this);};


        }
        exporter.close();

        if(visuals != null) visuals.close();

        // Close the simulation (PROBABLY NEVER REACHED)
        this.Close(logged);




    }


    // ------------------------------------------------------------------------------------------------------------
    // Manage and save output
    public void InitialiseCellLog(String cellCountLogFileName) {
        InitialiseCellLog(cellCountLogFileName, 1.);
    }

    public void InitialiseCellLog(String cellCountLogFileName, double frequency) {
        InitialiseCellLog(cellCountLogFileName,frequency,false);
    }

    public void InitialiseCellLog(String cellCountLogFileName, double frequency, Boolean profilingMode) {
        cellCountLogFile = new FileIO(cellCountLogFileName, "w");
        WriteLogFileHeader();
        this.cellCountLogFileName = cellCountLogFileName;
        this.logCellCountFrequency = frequency;
        if (profilingMode) {
            double[] tmpArr = GetModelState();
            int extraFields = extraSimulationInfoNames==null? 0: extraSimulationInfoNames.length;
            this.outputArr = new double[5][tmpArr.length+extraFields];
            // Initialise the logging array
            for (int i=0; i<outputArr.length; i++) {for (int j=0; j<outputArr[0].length; j++) {outputArr[i][j] = 0;}}
        }
    }

    private void WriteLogFileHeader() {
        cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells,DrugConcentration,rS,rR,mS,mR,dS,dR,dD,dt");
        if (extraSimulationInfoNames!=null) {
            cellCountLogFile.Write(",");
            cellCountLogFile.WriteDelimit(extraSimulationInfoNames, ",");
        }
        cellCountLogFile.Write("\n");

    }

    public void SetExtraSimulationInfo(String[] extraInfoNamesArr, double[] extraInfoArr) {
        this.extraSimulationInfoNames = extraInfoNamesArr;
        this.extraSimulationInfo = extraInfoArr;
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

    public void SaveTumourImage(int currTimeIdx) {

        if (imageFrequency > 0 && (currTimeIdx % (int) (imageFrequency/dt)) == 0) { //Save images at regualar intervals
            this.vis.ToPNG(mainDir.concat("/output/img/img_t_"+currTimeIdx*dt+".png"));
            resources.currV.ToPNG(mainDir.concat("/output/img/diff_img_t_"+currTimeIdx*dt+".png"));
            System.out.print("Saving image");
        } else if(imageFrequency<0 && currTimeIdx ==1 ){ // save image only in the beginning and end
            this.vis.ToPNG(mainDir.concat("/output/img/img_t_"+currTimeIdx*dt+".png"));
            resources.currV.ToPNG(mainDir.concat("/output/img/diff_img_t_"+currTimeIdx*dt+".png"));

        }else if(imageFrequency<0 && currTimeIdx ==tEnd ){ // save image at the end
            this.vis.ToPNG(mainDir.concat("/output/img/img_t_"+currTimeIdx*dt+".png"));
            resources.currV.ToPNG(mainDir.concat("/output/img/diff_img_t_"+currTimeIdx*dt+".png"));
        }
    }


    public void SaveModelState(String stateFileName) {
        // Can't have active pointers when saving the model. So, close everything here.
        if (cellCountLogFile!=null) {
            cellCountLogFile.Close();
            cellCountLogFile = null;
        }
        if (vis!=null) {vis = null;}
        SaveState(this,stateFileName);
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