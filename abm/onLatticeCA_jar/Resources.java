package onLatticeCA_jar;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.UIGrid;
import static HAL.Util.HeatMapRGB;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * This class managers blood and resource diffusion
 * Loads up user-fed vessel locations
 * supervises blood diffusion
 */
public class Resources {
    OnLatticeGrid onlatticeGrid;
     UIGrid currV;
     PDEGrid2D pdegrid2d;
     public double maximum_delta;
     double non_dim_diff_rate;
     int[][] vessels;
     int[] vessel_indices;


    public Resources(OnLatticeGrid G, int xdim, int ydim) {
        this.onlatticeGrid=G;
        currV = new UIGrid(G.xDim, G.yDim,  G.scaleFactor, G.visualiseA);
        non_dim_diff_rate=(G.Params.Resource_diff_rate*G.dt)/Math.pow(Params.spaceStep,2);
        pdegrid2d = new PDEGrid2D(xdim, ydim,true,true);

    }
    public Resources(){

    }
    public void initialize_vessels(String path2Vessels){
        List<int[]> rowsList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path2Vessels))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split("\\s+");  // Split by spaces or tabs
                int[] row = new int[values.length];    // Create a row array

                // Convert each string value to an integer
                for (int i = 0; i < values.length; i++) {
                    row[i] = Integer.parseInt(values[i]);
                }

                rowsList.add(row);  // Add the row to the list
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Convert the list to a 2D array
        vessels = rowsList.toArray(new int[rowsList.size()][]);
    }

    public  void Drawvessels(){
        for (int x = 0; x < onlatticeGrid.xDim; x++) {
            for(int y=0;y<onlatticeGrid.yDim; y++) {
                currV.SetPix(x, y,  HeatMapRGB(pdegrid2d.Get(x, y)/Params.vssl_bdary_value));// Scale for visual purposes
            }
        }
    }


    public void vessels2Indices(){ //Converts the vessel coordinates to position idexes
        int[] vesselIndex  = new int[vessels.length];
        for (int i = 0; i < vessels.length; i++) {
            vesselIndex[i] = onlatticeGrid.I(vessels[i][0], vessels[i][1]); //convert vessel coordinates to index
        }
        vessel_indices=vesselIndex;
    }

    public double setDirechletCond() {
        // advisable to limit Courant number to 0.2/
        double timescalar = Params.Courant_number/non_dim_diff_rate;
        pdegrid2d.DiffusionADI(non_dim_diff_rate*timescalar);
        for(Cell c:onlatticeGrid){
            c.Consumption(timescalar);
        }
        for (int[] vessel : vessels) {
            pdegrid2d.Set(vessel[0], vessel[1], Params.vssl_bdary_value);
        }
         maximum_delta = pdegrid2d.MaxDeltaScaled(0.);
        pdegrid2d.Update();
        return timescalar;
    }





}
