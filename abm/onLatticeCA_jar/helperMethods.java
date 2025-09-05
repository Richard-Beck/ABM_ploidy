package onLatticeCA_jar;

import HAL.GridsAndAgents.PDEGrid2D;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class helperMethods {
    /**
     * This class contains miscellaneous methods for various uses eg. printing Matrices to the console,
     * finding relative complement of two int arrays etc.
     */
    // Resources resources=new Resources();

    //Ensure copy numbers are within a  bound e.g [1,4].
    public static boolean isAboveThreshold(int limit, int[] data) {
        for (int k = 0; k < data.length; k++) {
            if (data[k] >= limit || data[k] <= 0)
                return true;
        }
        return false;
    }

    public static List<List<Double>> readMatrixFromFile(String filePath) throws IOException {
        List<List<Double>> matrix = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue; // Skip empty lines
                }

                List<Double> row = new ArrayList<>();
                String[] values = line.split("\\s+");
                for (String value : values) {
                    try {
                        row.add(Double.parseDouble(value));
                    } catch (NumberFormatException ex) {
                        System.err.println("Skipping invalid numeric value: " + value);
                    }
                }
                if (!row.isEmpty()) {
                    matrix.add(row);
                }
            }
        }

        if (matrix.isEmpty()) {
            System.err.println("Warning: No valid data found in file: " + filePath);
        }

        return matrix;
    }

    public boolean IsMember(int index, int[] array){
        helperMethods helper=new  helperMethods();
        boolean test=false;
        for (int element : array) {
            if (element == index) {
                test = true;

            }
        }
        return test;
    }

    public static double getTotalO2ConcInTheGrid(PDEGrid2D grid) { // Gets the total oxygen in the grid
        int xDim=grid.xDim;
        int yDim=grid.yDim;
        double totalO2Conc=0;
        for(int i=0;i<xDim;i++){
            for(int j=0;j<yDim;j++){
                totalO2Conc += grid.Get(i,j);
            }
        }

        return totalO2Conc;
    }


}


