package onLatticeCA_jar;

import HAL.Util;

import java.io.*;
import java.util.*;
import java.lang.reflect.Field;

public class Params {

    /**
     * This class managers all parameter values and reads input data from file
     */
    OnLatticeGrid G;

    String path2Vessels;

    //Resources resources=new Resources();
    final static String fileOut = "out/PopOut.csv";
    final static int BLACK= Util.RGB(0,0,0);
    final static int BLUE= Util.RGB(0,0,1);
    final static int RED= Util.RGB(1,0,0);
    final static int GREEN= Util.RGB(0,1,0);
    final static int YELLOW= Util.RGB(1,1,0);
    final static int CYAN= Util.RGB(0,1,1);


    int InitialPopulationType; //Selects the type of initialization, either seed from real data 0, or randomly seed diploi cells
    //Grid Parameters
    static public int xDim;// x-dimension of domain (in lattice sites)
    static public int yDim; // y-dimension of domain (in lattice sites)
    static double spaceStep=0.02;
    static int seed;



    //Cell rates

    static double initialTumorSize;//=3000;
    static int[][] inputKaryotype;
    static double initial_oxygen_conc;//=2.973e-10;
    static double divisionRate_S = 0.057;
    static double divisionRate_R = 0.057;
    static double movementRate_S = 0;
    /**  Uncomment the following to set sensitive-resistant cell specific rates
     static double movementRate_R = 0;
     static double deathRate_S =0.0027; // Natural death rate of sensitive cells in d^-1. Mean life span of a cell will be 1/k
     static double deathRate_R = 0.0027;// Natural death rate of resistant cells in d^-1.
     */
    static double drugKillProportion = 0;//0.75; // Drug induced death rate in d^-1.
    static double misRateWeight; //Weight for computing missegregation rate of cells
    static double natural_mis_seg_rate;
    static double cancer_cell_mis_seg_rate;
    static double normal_cell_mis_seg_rate;
    //static double natural_mis_seg_rate_cancer=0.0075;
    //static double[] misSegRateValues={0.005};
    //static int[] tEndValues={300};

    static double  natural_death_rate;//=0.031;
    static double  natural_growth_rate;//=0.063;
    static double cancer_cell_growth_rate;

    static int runGUI = 1;



    //Oxygen Diffusion parameters
    static double Resource_diff_rate;//Math.pow(10,-5);;//Math.pow(10,-5);// Math.pow(10,-5) 0.864 cm2 per day
    static double Resource_flux; //=100;//9.98e-12;//.25;//8.622*Math.pow(10,-7);

    //Following array for running replicates
    //    static int[] flow_levels={3250, 3000,2750,2500,2250,2000,1750,1500,1250,1000,750,500,250};
    //    static int[] initial_levels={100,80,50,20,10,0};
    static int[] oxygenAtVessels={4000};
    static int[] InitialSizes={20000,10000,5000};
    static double[] consumptionrates={100.0,30.0,20.0, 15.0, 10.0, 5.0,1.0};
    static double vssl_bdary_value;
    static double IniO2ConcScaling=Math.pow(10,-3)/2;//MTB .0001 to account for 10^4 conversion, and
    static double diff_vis_scale;//Math.pow(10,1);//  This scales the small oxygen values for visualization purposes
    static double normal_oxygen_consumption_rate;
    // MTB 10^-3 to convert to mg/cm^3
    static double cancer_oxygen_consumption_rate;

    static boolean Resource_on;
    static int Resource_start=1;
    static int Resource_stop;
    static int diff_time_step; //How many times diffusion and consumption is run per unit cell time
    static int con_time_step;

    static double diffusion_dt;

    static double diffusion_tol;

    static int writeFrequency;
    static int imageFrequency;



    //Other rates
    static double dt;// = 1; // Time step in d
    static int tEnd;// simulation end time
    static double consumption_MM_constant; //Consumption Michaelis Menten constant
    static double growth_michaelis_menten_constant;
    static double death_michaelis_menten_constant;

    //Data writing
    static boolean writeO2Data=false; //write oxygen concentration data per grid cell and total concentration in the tissue at specified  time points
    static boolean writeKaryotypeData=true; //Turn on/off when to write karyotype data per grid cell at specified  time points
    static boolean writeCellTypeData=false;   //Turn on/off when to write regional data. Data on Cells in a specified local disks with a given radius
    static boolean writeDistanceDensityData=false; //Turn on/off when to write regional data. Spot-wise Data on dist to vessels and density

    int Radius=3;
    static int Spot=200;  // in the case of multiple spots, Selects a file associated with a certain random set of spots




    static int nvessels; //number of blood vessels
    static int[][] bloodVessels;




    int[][] karyoLocationMatrix;
    //holds user input data (n by 24) comprising cell location with associated karyotype




    public Params(){

    }

    //Load parameters
    public void parameter(String path2Parameters) {
        File parameters = new File(OnLatticeGrid.mainDir.concat(path2Parameters));
        FileReader file;

        {
            try {
                file = new FileReader(OnLatticeGrid.mainDir.concat(path2Parameters));
                BufferedReader bf = new BufferedReader(file);
                String st = bf.readLine();
                int sum = 0;
                while ((st = bf.readLine()) != null) {
                    StringTokenizer stn = new StringTokenizer(st);
                    String rn = stn.nextToken();
                    String name = stn.nextToken();
                    int phy = Integer.parseInt(stn.nextToken());
                    System.out.println("Parameter " + rn);
                }
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }


    double[][] GrowthDeathRates;

    {
        try {
            GrowthDeathRates = readRates();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    int[][] centers;//  Centers of chosen spots for analysis after simulation
    {
        try {
            centers = testSpots();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    //

    public void load(String filename) {
        Map<String, String> params = new HashMap<>();

        // Read the file line by line
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\\s+"); // Split by spaces or tabs
                if (parts.length == 2) {
                    params.put(parts[0], parts[1]);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Update the fields using reflection
        for (Map.Entry<String, String> entry : params.entrySet()) {
            try {
                Field field = this.getClass().getDeclaredField(entry.getKey());
                field.setAccessible(true);

                // Check the field type and parse accordingly
                if (field.getType() == int.class) {
                    field.setInt(this, Integer.parseInt(entry.getValue()));
                } else if (field.getType() == double.class) {
                    field.setDouble(this, Double.parseDouble(entry.getValue()));
                } else if (field.getType() == boolean.class) {
                    field.setBoolean(this, Boolean.parseBoolean(entry.getValue()));
                } else if (field.getType() == String.class) {
                    field.set(this, entry.getValue());
                }
            } catch (NoSuchFieldException | IllegalAccessException e) {
                System.out.println("Warning: No matching field found for parameter: " + entry.getKey());
            }
        }
    }


  public int[][] readVessel() throws IOException {
      List<List<Double>> ves = helperMethods.readMatrixFromFile(G.mainDir.concat("/input/vessels.txt"));//Reads vessel  data as list
       bloodVessels = ves.stream().map(x -> x.stream().mapToInt(Double::intValue).toArray()).toArray(int[][]::new); //converts vessel data to matrix

      if(bloodVessels.length==1 && bloodVessels[0].length!=0){

          this.nvessels =1;
      } else{
          this.nvessels=bloodVessels.length;
      }


      //helperMethods.printMatrix(bloodVessels);
      //System.out.println("nvessels "+nvessels);


      List<List<Double>> kar = null;//Reads vessel  data as list
      try {
          kar = helperMethods.readMatrixFromFile(G.mainDir.concat("/input/ctrl_karyotype.txt"));
      } catch (IOException e) {
          throw new RuntimeException(e);
      }
      karyoLocationMatrix = kar.stream().map(x -> x.stream().mapToInt(Double::intValue).toArray()).toArray(int[][]::new); //converts location data to matrix
      this.inputKaryotype=karyoLocationMatrix;
   return  karyoLocationMatrix;
  }



    public int[][]  testSpots() throws IOException {

    List<List<Double>> spots = helperMethods.readMatrixFromFile(G.mainDir.concat("/input/testSpots/testSpot"+Spot+".txt"));//Reads coordinates for spots as list
    int[][] centers  = spots.stream().map(x -> x.stream().mapToInt(Double::intValue).toArray()).toArray(int[][]::new); //converts spot data to matrix
    return centers;

}

    public double[][] readRates() throws IOException {

        List<List<Double>> rat = helperMethods.readMatrixFromFile(G.mainDir.concat("/input/GrowthDeathRates.txt"));

        double[][] rates =rat.stream().map(x -> x.stream().mapToDouble(Double::doubleValue).toArray()).toArray(double[][]::new);
        return rates;
    }


}
