package onLatticeCA_jar;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.UIGrid;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.util.List;
import java.util.Arrays;
import java.util.Set;


import static HAL.Util.ArrayMax;
import static HAL.Util.HeatMapRGB;


/**
 * This class managers blood and resource diffusion
 * Loads up user-fed vessel locations
 * supervises blood diffusion
 */
public class Resources {
    helperMethods helper=new helperMethods();
      OnLatticeGrid onlatticeGrid;

     UIGrid currV;

      //Rand rn;
     PDEGrid2D pdegrid2d;

     FileIO vesselReader = null;
    String inputDir = "/Users/4477116/Documents/hal_wrapper/input"; // Directory  to read data
    //public static UIGrid vis2;
    public  int numberOfVessels;
    public double boundayValue;
    public double maximum_delta;
    public double[] maximum_delta_cus;



    /**

    static List<List<Integer>> vessselsData;//Reads vessel  data as list from txt

    {
        try {
            vessselsData = helper.readMatrixFromFile(onlatticeGrid.mainDir.concat("/input/vessels.txt"));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    static int[][] bloodVessels  = vessselsData.stream().map(x -> x.stream().mapToInt(Integer::intValue).toArray()).toArray(int[][]::new); //con
     */




    public Resources(OnLatticeGrid onlatticeGrid, int xdim, int ydim) {
        this.onlatticeGrid=onlatticeGrid;
        currV = new UIGrid(onlatticeGrid.xDim, onlatticeGrid.yDim,  onlatticeGrid.scaleFactor, onlatticeGrid.visualiseA);
        //public MelanomaModel(int x, int y, UIGrid vis, String outputFileName) {
       // super();
        pdegrid2d = new PDEGrid2D(xdim, ydim,true,true);

    }
    public Resources(){

    }
    static double non_dim_diff_rate=(Params.Resource_diff_rate*0.11)/Math.pow(Params.spaceStep,2); //Calculating the nondimensionalized diffusion coefficient
    //Setup vessels
    public int[][] bloodVessels;
    public int actual_nvessels;
        //static int nvessels=bloodVessels.length;






     static int[][] vessels=new int[Params.nvessels][2];//=new int[Params.nvessels][2];



    public void InitialO2Concentrations() throws Exception {
        //Use to read specific user initial oxygen conccentration to initialize simulation
        List<List<Double>> ves = helperMethods.readMatrixFromFile(onlatticeGrid.mainDir.concat("/input/oxygen.txt"));//Reads vessel  data as list
        int[][] O2  = ves.stream().map(x -> x.stream().mapToInt(Double::intValue).toArray()).toArray(int[][]::new); //converts vessel data to matrix

        double[][] InitialO2Conc=new double[O2.length][O2.length];

        for(int i=0;i<O2.length;i++){
            for(int j=0;j<O2.length;j++){
                InitialO2Conc[i][j]=Math.abs(O2[i][j]);
                //InitialO2Conc[i][j]=15*Math.pow(10,-9);
            }


        }


        // set oxygen levels in the grid
        for(int i=0;i<Params.xDim;i++){
            for(int j=0;j<Params.yDim;j++){

                //pdegrid2d.Set(i,j,InitialO2Conc[i][j]); //In the even that location specific O2 concetration values are used
                pdegrid2d.Set(i,j,Params.initial_oxygen_conc); //Constant initial concentrations at all locations
            }

        }

        pdegrid2d.Update();

    }

    public int[] InitVessels() throws IOException {
        int[] vesselIndex  = new int[Params.nvessels]; // Contains corresponding  position indexes for vessel coordinates (x,y)

         for(int i=0;i< Params.nvessels;i++){
          // helperMethods.printMatrix(vessels);
            vessels[i][0]=Params.bloodVessels[i][0];
            vessels[i][1]=Params.bloodVessels[i][1];
             vesselIndex[i]=onlatticeGrid.I(Params.bloodVessels[i][0],Params.bloodVessels[i][1]);
        }


        /**
        Uncomment to set vessels randomly

        rn=new Rand();
        for(int i=0; i<Params.nvessels; i++){
            vessels[i][0] = rn.Int( onlatticeGrid.xDim);
            vessels[i][1] = rn.Int( onlatticeGrid.yDim);
        }
        helperMethods.printMatrix(vessels);
         */
        return vesselIndex;

    }

//@TODO covert and save vessel index once
    public int[] vesselIndex(){ //Converts the vessel coordinates to position idexes
        int[] vesselIndex  = new int[vessels.length];
        for (int i = 0; i < vessels.length; i++) {
            vesselIndex[i] = onlatticeGrid.I(vessels[i][0], vessels[i][1]); //convert vessel coordinates to index
        }

        return  vesselIndex;
    }

    public  void Drawvessels(){
        for (int x = 0; x < onlatticeGrid.xDim; x++) {
            for(int y=0;y<onlatticeGrid.yDim; y++) {
                currV.SetPix(x, y,  HeatMapRGB(pdegrid2d.Get(x, y)/Params.vssl_bdary_value));// Scale for visual purposes
            }
        }
    }

    public double setDirechletCond() {
        // advisable to limit Courant number to 0.2/
        double timescalar = 0.2/non_dim_diff_rate;
        pdegrid2d.DiffusionADI(non_dim_diff_rate*timescalar);
        for(Cell c:onlatticeGrid){
            c.Consumption(timescalar);
        }
        for (int i = 0; i < Params.nvessels; i++) {
            pdegrid2d.Set(vessels[i][0], vessels[i][1],Params.vssl_bdary_value);
        }
        double[] deltas = pdegrid2d.GetDeltas();
        maximum_delta = 0;
        for(int i = 0; i< deltas.length; i++){
            maximum_delta=Math.max(maximum_delta,deltas[i]/pdegrid2d.Get(i));
        }
        //maximum_delta = pdegrid2d.MaxDelta();
        pdegrid2d.Update();
        return timescalar;
    }

    public void DiffAll(boolean Resource_on) { // This method is the flux version of setDirichletCond above. instead of a fix value at vessel site,
        //this method add a flux at set time points

        if (Resource_on) {

            //Resources  will enter through boundaries
            for (int i = 0; i < Params.nvessels; i++) {
                //double inflow=pdegrid2d.Get(vessels[i][0],vessels[i][1]);
                double diff=Params.vssl_bdary_value-pdegrid2d.Get(vessels[i][0], vessels[i][1]);
                pdegrid2d.Add(vessels[i][0], vessels[i][1],(Params.Resource_flux));
                //pdegrid2d.Add(vessels[i][0], vessels[i][1], (Params.Resource_flux));
  /**
                if(diff>Params.Resource_flow_value) {
                    pdegrid2d.Add(vessels[i][0], vessels[i][1], (Params.Resource_flow_value));
                }else{
                    pdegrid2d.Add(vessels[i][0], vessels[i][1], diff);
                }
  */
            }

        } else {

            for (int i = 0; i < Params.nvessels; i++) {
                pdegrid2d.Add(vessels[i][0], vessels[i][1], 0);
            }

    }
        pdegrid2d.DiffusionADI(non_dim_diff_rate);
        /**
        for(Cell c:onlatticeGrid){
            //c.Consumption(c); //This is commented out and the consumption method below is used

        }
        */
        pdegrid2d.Update();

    }




}