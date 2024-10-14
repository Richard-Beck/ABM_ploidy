package onLatticeCA_jar;


import HAL.Rand;

import java.util.Arrays;

import java.util.stream.IntStream;

import static HAL.Util.*;

/**
 * This class handles initialization of the cell population.
 * Loads user data on initial population comprising cell locations (x,y) and karyotype (22-element vec)
 * To use  initial population data, set "InitialPopLocation" to zero in the main class
 * To randomize initial population and assign same karyotype (diploid) set the value of "InitialPopLocation" between 1 and 5
 */


import static HAL.Util.HeatMapRGB;
public class Initialization {
    OnLatticeGrid G;
    helperMethods helper = new helperMethods();
    Rand rand;
    Params params;

public Initialization(OnLatticeGrid G){
    this.G=G;

}

    public int[] getIniLocationIndex() { //user input

        int[] vesselIndex = G.resources.vesselIndex();
        int[][] Locations = Params.inputKaryotype;//this.getInputKaryotype();
        //int[][] initialPopKaryotype =new int[Locations.length][22]; //Hold karyotypes information
        int[] initialPopLocIndex = new int[Locations.length]; //holds Location indexes

        for (int i = 0; i < Locations.length; i++) {
            initialPopLocIndex[i] = G.I(Locations[i][0], Locations[i][1]);

        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialPopLocIndex = helper.relativeComplement(initialPopLocIndex, vesselIndex);

        if(viableInitialPopLocIndex.length<Locations.length){
            //If cell's location happens to be the same as a vesel location, cell is assigned a different location
           int difference=Locations.length-viableInitialPopLocIndex.length;

            int[] AllLocation = IntStream.rangeClosed(1, (G.xDim*G.yDim)).toArray();
            int[]unusableLocation=helper.getUnion(viableInitialPopLocIndex,vesselIndex);
            int[] usableLocation=helper.relativeComplement(AllLocation,unusableLocation);
            int[]randomIndex=helperMethods.selectRandomSubset(usableLocation,difference);
            System.out.println("Random "+Arrays.toString(randomIndex));

            viableInitialPopLocIndex=helper.getUnion( viableInitialPopLocIndex,randomIndex);

        }
        System.out.println("All locations "+initialPopLocIndex.length);
        System.out.println("Viable locations left: " + viableInitialPopLocIndex.length);
        System.out.println("initialPopLocIndex "+Arrays.toString(initialPopLocIndex));

        System.out.println("Viable locations "+Arrays.toString(viableInitialPopLocIndex));


        return viableInitialPopLocIndex;
    }


    public int[][] getIniPopKaryotype() { // sets initial cell locations from user input
        int[] vesselIndex = G.resources.vesselIndex();
        int[][] Locations = Params.inputKaryotype;//this.getInputKaryotype();
        int[][] initialPopKaryotype = new int[Locations.length][22]; //Hold karyotypes information
        for (int i = 0; i < Locations.length; i++) {
            for (int j = 2; j < Locations[1].length; j++) {
                initialPopKaryotype[i][j - 2] = Locations[i][j];
            }
        }

        return initialPopKaryotype;

    }

    public int[] InitialCenterCirc() { // Gets viable grid cell locations in the center of the grid as disk for seeding
        int[] vesselIndex = G.resources.vesselIndex();
        double Intended_radius=Math.ceil((double)Math.sqrt(Params.initialTumorSize/Math.PI));
        int gridRadius=Math.min(Params.xDim/2,Params.yDim/2);
        double correction=Math.ceil(Math.max((Intended_radius-gridRadius),0));
        int Radius=(int)(Intended_radius+correction);
        // Function to select values of initialLocArr that belong to the center of the grid
        int[] tumorNeighborhoodCenter = CircleHood(true, Radius); //indexes for centered circular initialization
        int hoodSize = G.MapHood(tumorNeighborhoodCenter, G.xDim / 2, G.yDim / 2);
        int[] initialLocArro = tumorNeighborhoodCenter;

        /** initialLocArro contains all indexes followed by their respective  (x,y) coordinates.
        We are only interested in the indexes for now, so we delete x-y coordinates
        */
        int n = (initialLocArro.length / 3);
        int[] viableInitialLocations = new int[n];
        for (int i = 0; i < n; i++) {
            viableInitialLocations[i] = initialLocArro[i];
        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(viableInitialLocations, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;
    }

    public int[] InitialTopRightCirc() { //Seeds cells in the top right conner in the form of a disk
        int[] vesselIndex = G.resources.vesselIndex();
        double Intended_radius=Math.ceil((double)Math.sqrt(Params.initialTumorSize/Math.PI));
        int gridRadius=Math.min(Params.xDim/2,Params.yDim/2);
        double correction=Math.ceil(Math.max((Intended_radius-gridRadius),0));
        int Radius=(int)(Intended_radius+correction);
        int[] tumorNeighborhoodTopRight = CircleHood(true,Radius); //indexes for top right circular initialization
        int hoodSizeTopRight = G.MapHood(tumorNeighborhoodTopRight, (G.xDim - 25), (G.yDim - 25));
        int[] initialLocArro = tumorNeighborhoodTopRight;

        // initialLocArro containing indexes followed by their (x,y) coordinates.
        //We are only interested in the indexes for now, so we delete x-y coordinates
        int n = (initialLocArro.length / 3);
        int[] viableInitialLocations = new int[n];
        for (int i = 0; i < n; i++) {
            viableInitialLocations[i] = initialLocArro[i];
        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(viableInitialLocations, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;

    }


    public int[] InitialTopLeftCirc() { //Seeds cells in the top left conner in the form of a disk
        int[] vesselIndex = G.resources.vesselIndex();
        int[] tumorNeighborhoodTopLeft = CircleHood(true, Math.sqrt(Params.initialTumorSize/Math.PI)+5); //indexes for top left circular initialization
        int hoodSizeTopLeft = G.MapHood(tumorNeighborhoodTopLeft, 25, (G.yDim - 25));

        int[] initialLocArro = tumorNeighborhoodTopLeft;

        // initialLocArro containing indexes followed by their (x,y) coordinates.
        //We are only interested in the indexes for now, so we delete x-y coordinates
        int n = (initialLocArro.length / 3);
        int[] viableInitialLocations = new int[n];
        for (int i = 0; i < n; i++) {
            viableInitialLocations[i] = initialLocArro[i];
        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(viableInitialLocations, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;

    }

    public int[] InitialBottomRightCirc() {
        int[] vesselIndex = G.resources.vesselIndex();
        int[] tumorNeighborhoodBottomRight = CircleHood(true, Math.sqrt(Params.initialTumorSize/Math.PI)+5); //indexes for bottom right circular initialization
        int hoodSizeBottomRight = G.MapHood(tumorNeighborhoodBottomRight, (G.xDim - 25), 25);


        int[] initialLocArro = tumorNeighborhoodBottomRight;

        // initialLocArro containing indexes followed by their (x,y) coordinates.
        //We are only interested in the indexes for now, so we delete x-y coordinates
        int n = (initialLocArro.length / 3);
        int[] viableInitialLocations = new int[n];
        for (int i = 0; i < n; i++) {
            viableInitialLocations[i] = initialLocArro[i];
        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(viableInitialLocations, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;

    }


    public int[] InitialBottomLeftCirc() {
        int[] vesselIndex = G.resources.vesselIndex();
        int[] tumorNeighborhoodBottomLeft = CircleHood(true, 21); //indexes for bottom right circular initialization
        int hoodSizeBottomLeft = G.MapHood(tumorNeighborhoodBottomLeft, 25, 25);


        int[] initialLocArro = tumorNeighborhoodBottomLeft;

        // initialLocArro containing indexes followed by their (x,y) coordinates.
        //We are only interested in the indexes for now, so we delete x-y coordinates
        int n = (initialLocArro.length / 3);
        int[] viableInitialLocations = new int[n];
        for (int i = 0; i < n; i++) {
            viableInitialLocations[i] = initialLocArro[i];
        }
        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(viableInitialLocations, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;

    }

    public int[] InitialSquareCenter() {
        int[] vesselIndex = G.resources.vesselIndex();
        int[][] initialLocMatrix = new int[G.xDim/2][G.yDim/2];
        for (int i = (G.xDim / 2 - 25); i < (G.xDim / 2 + 25); i++) {
            for (int j = (G.yDim / 2 - 25); j < (G.yDim / 2 + 25); j++) {
                int a = i - (G.xDim / 2 - 25);
                int b = j - (G.yDim / 2 - 25);
                initialLocMatrix[i - (G.xDim / 2 - 25)][j - (G.yDim / 2 - 25)] = G.I(i, j); //I(i,j) Converts the (i,j) coordinate to position index
                int c = G.I(i, j);
                int d = c + a;

            }


        }

        int[] initialLocArro  = helper.Matrix2Vector(initialLocMatrix); //indexes for centered rectangular initialization
        int[] viableInitialLocArro = helper.relativeComplement(initialLocArro, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;
    }


    public int[] InitialRandom() {
        int[] vesselIndex = G.resources.vesselIndex();
        int[] initialLocArr = new int[G.xDim * G.yDim];
        // Generate a list of random grid positions
        for (int i = 0; i < G.xDim * G.yDim; i++) {
            initialLocArr[i] = i;
        }
        helper.shuffle(initialLocArr);


        // Check if intended cell location does not overlap with a vessel location,
        int[] viableInitialLocArro = helper.relativeComplement(initialLocArr, vesselIndex);
        helper.shuffle(viableInitialLocArro);// Shuffle them up
        return viableInitialLocArro;

    }



    public int[] setLocation() {

        /**
         *Selecting  possible cell initialization location
         */
        int[] InitialUserInputLoc=this.getIniLocationIndex(); //Specific to user input data
        int[][] initialPopKaryotype =this.getIniPopKaryotype();  //Specific to user input data
        int[] InitialCenterCirc=InitialCenterCirc();
        int[] InitialTopRightCirc=InitialTopRightCirc();
        int[] InitialTopLeftCirc=InitialTopLeftCirc();
        int[] InitialBottomRightCirc=InitialBottomRightCirc();
        int[] InitialBottomLeftCirc=InitialBottomLeftCirc();
        int[] InitialSquareCenter=InitialSquareCenter();
        int[] InitialRandom=InitialRandom();

        int[] initialLocArro=InitialRandom;
        int [][] InitialKaryotype=G.baseKaryotype;
         //Using user input value for initialLocation, we decide on which region initial cell population will be placed
         if (G.InitialPopulationType == 0) {
             initialLocArro = InitialUserInputLoc;
         } else if (G.InitialPopulationType== 1) {
             initialLocArro = InitialCenterCirc;
         } else if (G.InitialPopulationType == 2) {
             initialLocArro = InitialTopRightCirc;
         } else if (G.InitialPopulationType == 3) {
            initialLocArro = InitialTopLeftCirc;
         } else if (G.InitialPopulationType == 4) {
              initialLocArro = InitialBottomRightCirc;
         } else if (G.InitialPopulationType == 5) {
              initialLocArro = InitialBottomLeftCirc;
         } else if (G.InitialPopulationType == 6) {
             initialLocArro = InitialSquareCenter;
         }else if (G.InitialPopulationType == 7) {

             initialLocArro = InitialRandom;
         }

        return initialLocArro;

    }

}
