package onLatticeCA_jar;

import HAL.GridsAndAgents.PDEGrid2D;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Random;
import java.util.Arrays;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.List;
import java.text.SimpleDateFormat;
import java.util.Date;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.stream.IntStream;

import static HAL.Util.MooreHood;


public class helperMethods {
    /**
     * This class contains miscellaneous methods for various uses eg. printing Matrices to the console,
     * finding relative complement of two int arrays etc.
     */
    // Resources resources=new Resources();
    public static Set<Integer> arrayUnion(int[] a, int[] b) {
        int n = a.length;
        int m = b.length;
        // find min of n and m
        int min = (n < m) ? n : m;

        // set container
        Set<Integer> set = new HashSet<>();

        // add elements from both the arrays for
        // index from 0 to min(n, m)-1
        for (int i = 0; i < min; i++) {
            set.add(a[i]);
            set.add(b[i]);
        }

        // add remaining elements to the set from the other
        // array (having greater length)
        // note that only one of the loops will execute
        if (n > m) {
            for (int i = m; i < n; i++) {
                set.add(a[i]);
            }
        } else if (n < m) {
            for (int i = n; i < m; i++) {
                set.add(b[i]);
            }
        }

        // driver code to print the output


        //return set();

        // return  helperMethods.getUnion();

        return set;
    }

    // Driver Code
    public static int[] getUnion(int[] a, int[] b) {
        Integer[] temp = helperMethods.arrayUnion(a, b).toArray(new Integer[0]);
        //System.out.print("temp: "+ Arrays.toString(temp));
        int[] output = new int[temp.length];
        for (int i = 0; i < temp.length; i++) {
            output[i] = temp[i];
        }
        // System.out.print("Output: "+Arrays.toString(output));
        return output;
    }


    //Ensure copy numbers are within a  bound e.g [1,4].
    public static boolean isAboveThreshold(int limit, int[] data) {
        for (int k = 0; k < data.length; k++) {
            if (data[k] >= limit || data[k] <= 0)
                return true;
        }
        return false;
    }

    //To visualize 2D matrices
    public static void printMatrix(int matrix[][]) {
        for (int[] row : matrix) {
            // traverses through number of rows
            for (int element : row) {
                // 'element' has current element of row index
                System.out.print(element + "\t");
            }
            System.out.println();
        }
    }

    public static void printDoubleMatrix(double matrix[][]) {
        for (double[] row : matrix) {
            // traverses through number of rows
            for (double element : row) {
                // 'element' has current element of row index
                System.out.print(element + "\t");
            }
            System.out.println();
        }
    }

    public int[] Matrix2Vector(int[][] input) { // converts a matrix to a vector by row
        int[] out = new int[input.length * input[0].length];
        for (int i = 0; i < input.length; i++) {
            for (int j = 0; j < input[i].length; j++) {
                out[i + (j * input.length)] = input[i][j];
            }
        }
        return out;
    }


    public static Object Matrix2Vector(Object input) {
        int numRows = 0;
        int numCols = 0;

        // Determine the number of rows and columns based on the input type
        if (input instanceof int[][]) {
            int[][] intMatrix = (int[][]) input;
            numRows = intMatrix.length;
            numCols = intMatrix[0].length;
        } else if (input instanceof double[][]) {
            double[][] doubleMatrix = (double[][]) input;
            numRows = doubleMatrix.length;
            numCols = doubleMatrix[0].length;
        } else {
            throw new IllegalArgumentException("Unsupported input type");
        }

        // Convert the matrix to a vector by row
        Object out;
        if (input instanceof int[][]) {
            int[] intVector = new int[numRows * numCols];
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    intVector[i * numCols + j] = ((int[][]) input)[i][j];
                }
            }
            out = intVector;
        } else {
            double[] doubleVector = new double[numRows * numCols];
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    doubleVector[i * numCols + j] = ((double[][]) input)[i][j];
                }
            }
            out = doubleVector;
        }
        return out;
    }


    //converts a 1D vector of numbers to a single integer
    public static int array2String(int[] Array) {
        String sb = "";
        for (int i = 0; i < Array.length; i++) {
            sb = sb + (Integer.toString(Array[i]));
        }
        return Integer.valueOf(sb);
    }

    public static List<List<Integer>> readIntMatrixFromFile(String filePath) throws IOException { //Read matrix data from txt file
        List<List<Integer>> matrix = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                List<Integer> row = new ArrayList<>();
                String[] values = line.trim().split("\\s+");
                for (String value : values) {
                    try {
                        row.add(Integer.parseInt(value));
                    } catch (NumberFormatException ex) {
                    }// handle your exception
                }
                matrix.add(row);
            }
        }

        return matrix;
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

       /** //Print matrix
        for (List<Double> row : matrix) {
            System.out.println(row);
        }
      */
        return matrix;
    }






    public static void readAndAppend(int seed, String filePath) {

        // Read from the file
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            System.out.println("Content seed file:");
            while ((line = reader.readLine()) != null) {
                //System.out.println(line);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        // Append to the file
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) { // true for append mode
            String dateTime = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date());
            String textToAppend = seed + " " + dateTime;
            writer.write(textToAppend);
            writer.newLine(); // Add a new line after the appended text
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static void shuffleArray(int[] array, Random random) {
        for (int i = array.length - 1; i > 0; i--) {
            // Generate a random index between 0 and i (inclusive)
            int j = random.nextInt(i + 1);

            // Swap elements at indices i and j
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }
    }

    public static List<List<Integer>> readDoubleMatrixFromFile(String filePath) throws IOException { //Read matrix data from txt file
        List<List<Integer>> matrix = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                List<Integer> row = new ArrayList<>();
                String[] values = line.trim().split("\\s+");
                for (String value : values) {
                    row.add(Integer.parseInt(value));
                }
                matrix.add(row);
            }
        }

        return matrix;
    }



        public static double [][] Reading2DArrayFromFile(String path) throws Exception {
            Scanner sc = new Scanner(new BufferedReader(new FileReader(path)));
            int rows = 100;
            int columns = 100;
            double [][] myArray = new double[rows][columns];
            while(sc.hasNextLine()) {
                for (int i=0; i<myArray.length; i++) {
                    String[] line = sc.nextLine().trim().split(" ");
                    for (int j=0; j<line.length; j++) {
                        myArray[i][j] = Integer.parseInt(line[j]);
                    }
                }
            }
            System.out.println(Arrays.deepToString(myArray));
            return myArray;
        }


    private static char[][] finalmatrix;




    // Java program to find all those
// elements of arr1[] that are not
// present in arr2[]


        // Driver code

    public static int[] relativeComplement(int[] arr1, int[] arr2) { //Finds all elements in arr1 but not in arr2
        ArrayList<Integer> interest = new ArrayList<>();
        for (int i = 0; i < arr1.length; i++) {
            int toCheckValue = arr1[i];
            boolean test = false;
             for (int element : arr2) {
                if (element == toCheckValue) {
                    test = true;
                    break;
                }

            }
            if (test!=true){
                interest.add(arr1[i]);
            }

        }

        int[] myComplement = new int[interest.size()];
        Arrays.setAll(myComplement, interest::get);

        //System.out.println("MyComplement "+Arrays.toString(myComplement)); //Just checking to see if it correct
        return myComplement;
    }


    static ArrayList<Integer> removeNegative(int[] arr) {
        ArrayList<Integer> newArr = new ArrayList<Integer>();

        for (int x : arr) {
            if (x >= 0) {
                newArr.add(x);
            }
        }

        for (int x : newArr) {
            System.out.print(x + " ");
        }
        return newArr;
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

    static int[] shuffle(int[] arr)
    {
        int n=arr.length;
        // Creating a object for Random class
        Random r = new Random();

        // Start from the last element and swap one by one. We don't
        // need to run for the first element that's why i > 0
        for (int i = n-1; i > 0; i--) {

            // Pick a random index from 0 to i
            int j = r.nextInt(i+1);

            // Swap arr[i] with the element at random index
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
        return arr;
    }

    static int countDistinctElements(int[] arr)
    { //Return tehe total number of distinct elements in an array
        int n=arr.length;
        // First sort the array so that
        // all occurrences become consecutive
        Arrays.sort(arr);
        //Count distinct elements first
        int res = 0;

        // Pick all elements one by one
        for (int i = 0; i < n; i++) {
            int j = 0;
            for (j = 0; j < i; j++)
                if (arr[i] == arr[j])
                    break;

            // If not printed earlier,
            // then print it
            if (i == j)
                res++;
        }

       return res;
    }
    static int[] DistinctElements(int[] arr) {
        //Returns unique elements in an array
        int[] uniqueElements = Arrays.stream(arr).distinct().toArray();

        return uniqueElements;
    }

    static int[] getRowI(int[][] arr,int row_num) { //return the ith row of an int matrix
        if (arr == null || arr.length <= row_num || row_num < 0) {
            throw new IllegalArgumentException("Invalid input or row number out of bounds");
        }
        int ncols = arr[1].length;
        int[] Row = new int[ncols];
        for (int i = 0; i < ncols; i++) {
            Row[i] = arr[row_num][i];
        }

        return Row;

    }
    static double[] getRowIDouble(double[][] arr,int row_num) { //return the ith row of an int matrix
        if (arr == null || arr.length <= row_num || row_num < 0) {
            throw new IllegalArgumentException("Invalid input or row number out of bounds");
        }
        int ncols = arr[1].length;
        double [] Row = new double[ncols];
        for (int i = 0; i < ncols; i++) {
            Row[i] = arr[row_num][i];
        }

        return Row;

    }

    static int[] getcolumnI(int[][] arr,int col_num) {//return the ith column of an int  matrix
        int nrows = arr.length;
        int[] col = new int[nrows];
        for (int i = 0; i < nrows; i++) {
            col[i] = arr[i][col_num];
        }

        return col;

    }
    static double[] getcolumnIDouble(double[][] arr,int columnIndex) { //return the ith column of a double matrix
        int nrows = arr.length;
        double[] col = new double[nrows];
        for (int i = 0; i < nrows; i++) {
            col[i] = arr[i][columnIndex];
        }

        return col;

    }


    static double[][] MatirceDifference(double[][] A, double[][] B){
        System.out.println("second dimension "+A[1].length);
        double[][] output=new double[A.length][A[1].length];
        for(int m=0;m<A.length;m++){
            for(int n=0;n<A[1].length;n++){
                output[m][n]=A[m][n]-B[m][n];
            }
        }
    return output;
    }


    public static double getMaxValue(double[][] numbers) {
        double maxValue = numbers[0][0];
        for (int j = 0; j < numbers.length; j++) {
            for (int i = 0; i < numbers[j].length; i++) {
                if (numbers[j][i] > maxValue) {
                    maxValue = numbers[j][i];
                }
            }
        }
        return maxValue;
    }


    public static double getMinValue(Object numbers) {
        if (numbers instanceof int[][]) {
            int[][] intNumbers = (int[][]) numbers;
            int minValue = intNumbers[0][0];
            for (int j = 0; j < intNumbers.length; j++) {
                for (int i = 0; i < intNumbers[j].length; i++) {
                    if (intNumbers[j][i] < minValue) {
                        minValue = intNumbers[j][i];
                    }
                }
            }
            return minValue;
        } else if (numbers instanceof double[][]) {
            double[][] doubleNumbers = (double[][]) numbers;
            double minValue = doubleNumbers[0][0];
            for (int j = 0; j < doubleNumbers.length; j++) {
                for (int i = 0; i < doubleNumbers[j].length; i++) {
                    if (doubleNumbers[j][i] < minValue) {
                        minValue = doubleNumbers[j][i];
                    }
                }
            }
            return minValue;
        } else {
            throw new IllegalArgumentException("Unsupported input type");
        }
    }


    public static double sumMatrixElements(double [][] arr){
        double sum = 0.0;
        for (int i=0;i<arr.length;i++) {
            for (int j=0;j<arr[1].length;j++) {
               sum+=arr[i][j];
                }
             }
        return sum;
        }

    public static int vectorSum(int[] arr) {

        // method for sum of elements in an array
         int sum=0;
         // Iterate through all elements and add them to sum
            for (int i = 0; i < arr.length; i++) {
                sum += arr[i];
            }
            return sum;
        }

    public static double vectorSumDouble(double[] arr) {   // Sums two 1D arrays together


        double sum=0;
        // Iterate through all elements and add them to sum
        for (int i = 0; i < arr.length; i++) {
            sum += arr[i];
        }
        return sum;
    }

    static int[] toIntArray(List<Integer> integerList) { //Converts and Integer list to int array
        int[] intArray = new int[integerList.size()];
        for (int i = 0; i < integerList.size(); i++) {
            intArray[i] = integerList.get(i);
        }
        return intArray;
    }


    public static double euclideanDistance(int [] arr1, int [] arr2) {
        if (arr1.length != arr2.length) {
            throw new IllegalArgumentException("Vectors must have the same dimensionality");
        }

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            double diff = arr1[i] - arr2[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    public static double euclideanDistanceDouble(double[] arr1, double[] arr2) { //Finds the Eucleadean distance between to double arrays
        if (arr1.length != arr2.length) {
            throw new IllegalArgumentException("Vectors must have the same dimensionality");
        }

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            double diff = arr1[i] - arr2[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    public static double manhattanDistance(int[] arr1, int[] arr2) { //Finds the Manhattan distance between to double arrays
        if (arr1.length != arr2.length) {
            throw new IllegalArgumentException("Vectors must have the same dimensionality");
        }

        int sum = 0;
        for (int i = 0; i < arr1.length; i++) {
            int diff = Math.abs(arr1[i] - arr2[i]);
            sum +=  diff;
        }
        return sum;
    }

    public static double Distance(Object array1, Object array2, String kind) { //Finds the Euclidean or Manhattan distance between two 1D arrays

        double distance = 0;
        if (array1 instanceof int[] && array2 instanceof int[]) {
            int[] arr1 = (int[]) array1;
            int[] arr2 = (int[]) array2;

            if (arr1.length != arr2.length) {
                throw new IllegalArgumentException("Vectors must have the same dimensionality");
            }

            if (kind == "manhattan") {

                for (int i = 0; i < arr1.length; i++) {
                    double diff = Math.abs(arr1[i] - arr2[i]);
                    distance += diff;
                }
            } else if (kind == "euclidean") {
                for (int i = 0; i < arr1.length; i++) {
                    double diff = arr1[i] - arr2[i];
                    distance += diff * diff;
                }
                distance = Math.sqrt(distance);
            }

        } else if (array1 instanceof double[] && array2 instanceof double[]) {
            double[] arr1 = (double[]) array1;
            double[] arr2 = (double[]) array2;

            if (arr1.length != arr2.length) {
                throw new IllegalArgumentException("Vectors must have the same dimensionality");
            }
            if (kind == "manhattan") {

                for (int i = 0; i < arr1.length; i++) {
                    double diff = Math.abs(arr1[i] - arr2[i]);
                    distance += diff;
                }
            } else if (kind == "euclidean") {
                for (int i = 0; i < arr1.length; i++) {
                    double diff = arr1[i] - arr2[i];
                    distance += diff * diff;
                }
                distance = Math.sqrt(distance);
            }



        }

        return distance;
    }


    public static double vectorSum(Object array) {

        // method for finding the sum of elements in an array
        double sum = 0;
        if (array instanceof int[]) {
            int[] arr = (int[]) array;

            // Iterate through all elements and add them to sum
            for (int i = 0; i < arr.length; i++) {
                sum += arr[i];
            }
        }else if (array instanceof double[]) {
            double[] arr = (double[]) array;

            // Iterate through all elements and add them to sum
            for (int i = 0; i < arr.length; i++) {
                sum += arr[i];
            }
        }
        return sum;
    }

    public static int[][] deleteZeroRows(int[][] matrix) { //Delete zero rows from a 2D array and return a new array with remaining rows
        // Count non-zero rows
        int count = 0;
        for (int[] row : matrix) {
            boolean isZeroRow = true;
            for (int num : row) {
                if (num != 0) {
                    isZeroRow = false;
                    break;
                }
            }
            if (!isZeroRow) {
                count++;
            }
        }

        // Create a new matrix without zero rows
        int[][] result = new int[count][matrix[0].length];
        int index = 0;
        for (int[] row : matrix) {
            boolean isZeroRow = true;
            for (int num : row) {
                if (num != 0) {
                    isZeroRow = false;
                    break;
                }
            }
            if (!isZeroRow) {
                result[index++] = row;
            }
        }

        return result;
    }

    public static double[] calculatePairwiseDistances(int [][] matrix1, int[][] matrix2,String type) { //Computes the pairwise distances betweem rows of a matrix
        int numRows1 = matrix1.length;
        int numRows2 = matrix2.length;
        double[] distances = new double[numRows1*numRows2];

        int index = 0;
        for (int i = 0; i < numRows1; i++) {
            for (int j = 0; j < numRows2; j++) {
                distances[index++] = Distance(matrix1[i], matrix2[j],type);
            }
        }

        return distances;
    }

    public static int[][] getSubmatrix(int[][] matrix) {
        // Check if the matrix has zero rows
        if (matrix.length == 0) {
            return new int[0][0]; // Return an empty matrix
        }

        int numRows = matrix.length;
        int numCols = matrix[0].length;

        // Create submatrix with columns 2 through end
        int[][] submatrix = new int[numRows][Math.max(0, numCols - 2)];
        for (int i = 0; i < numRows; i++) {
            for (int j = 2; j < numCols; j++) {
                submatrix[i][j - 2] = matrix[i][j];
            }
        }

        return submatrix;
    }



    public static Object columnBindVectors(Object vector1, Object vector2) { //column bind; the equivalent to "cbind" in R
        // Determine the number of rows (size of vectors)
        int numRows = 0;
        if (vector1 instanceof int[]) {
            numRows = ((int[]) vector1).length;
        } else if (vector1 instanceof double[]) {
            numRows = ((double[]) vector1).length;
        } else {
            throw new IllegalArgumentException("Unsupported vector type");
        }

        // Create a new 2D array with two columns
        Object result;
        if (vector1 instanceof int[]) {
            int[][] intResult = new int[numRows][2];
            for (int i = 0; i < numRows; i++) {
                intResult[i][0] = ((int[]) vector1)[i]; // First column from vector1
                intResult[i][1] = ((int[]) vector2)[i]; // Second column from vector2
            }
            result = intResult;
        } else {
            double[][] doubleResult = new double[numRows][2];
            for (int i = 0; i < numRows; i++) {
                doubleResult[i][0] = ((double[]) vector1)[i]; // First column from vector1
                doubleResult[i][1] = ((double[]) vector2)[i]; // Second column from vector2
            }
            result = doubleResult;
        }

        return result;
    }





    public static Object extractStrictlyUpperTriangle(Object matrix) {//Extracts a strictly upper triangular matrix from a 2D array
        if (matrix instanceof int[][]) {
            int[][] intMatrix = (int[][]) matrix;
            int numRows = intMatrix.length;
            int numCols = intMatrix[0].length;
            int size = numRows * (numRows - 1) / 2;
            int[] upperTriangle = new int[size];
            int index = 0;
            for (int i = 0; i < numRows; i++) {
                for (int j = i + 1; j < numCols; j++) {
                    upperTriangle[index++] = intMatrix[i][j];
                }
            }
            return upperTriangle;
        } else if (matrix instanceof double[][]) {
            double[][] doubleMatrix = (double[][]) matrix;
            int numRows = doubleMatrix.length;
            int numCols = doubleMatrix[0].length;
            int size = numRows * (numRows - 1) / 2;
            double[] upperTriangle = new double[size];
            int index = 0;
            for (int i = 0; i < numRows; i++) {
                for (int j = i + 1; j < numCols; j++) {
                    upperTriangle[index++] = doubleMatrix[i][j];
                }
            }
            return upperTriangle;
        } else {
            throw new IllegalArgumentException("Unsupported input type");
        }
    }


    public Object Transpose(Object matrix) { //Finds the transpose of a matrix
        if (matrix instanceof double[][]) {
            double [][]arr=(double[][]) matrix;
            int rows = arr.length;
            int cols = arr[0].length;
            double[][] transpose = new double[cols][rows];
            for (int i=0;i < rows; i++) {

                for (int j = 0; j < cols; j++) {

                    transpose[j][i] = arr[i][j];
                }
            }
            return transpose;
        }else if(matrix instanceof int[][]){
            int [][]arr=(int[][]) matrix;
            int rows = arr.length;
            int cols = arr[0].length;

            int[][] transpose = new int[cols][rows];
            for (int i=0;i < rows; i++) {
                for (int j=0;j<cols; j++) {
                    transpose[j][i] = arr[i][j];

                }
            }
            return transpose;
        }else {
            throw new IllegalArgumentException("check input");
        }
    }



    public static double[][] unlistAndStackMatrices(List<double[][]> matrices) { // stacks all matrices from a list of matrices to form one big matrix
        List<double[]> elements = new ArrayList<>();
        for (double[][] matrix : matrices) {
            for (double[] row : matrix) {
                if (row.length != 2) {
                    throw new IllegalArgumentException("Each row should have 2 elements");
                }
                elements.add(row);
            }
        }

        double[][] result = new double[elements.size()][2];
        for (int i = 0; i < elements.size(); i++) {
            result[i] = elements.get(i);
        }

        return result;
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



    public static int[] selectRandomSubset(int[] array, int k) {
        Random rand = new Random();
        List<Integer> indices = new ArrayList<>();

        // Generate k unique random indices
        while (indices.size() < k) {
            int index = rand.nextInt(array.length);
            if (!indices.contains(index)) {
                indices.add(index);
            }
        }

        // Create subset array using selected indices
        int[] subset = new int[k];
        for (int i = 0; i < k; i++) {
            subset[i] = array[indices.get(i)];
        }

        return subset;
    }


   //used to read user input data for cell locations with associated karyotype
    public int[][] getInputKaryotype() {
        //Load user input initial cell population with karyotype information
        List<List<Double>> ves = null;//Reads vessel  data as list
        try {
            ves = helperMethods.readMatrixFromFile(OnLatticeGrid.mainDir.concat("/input/TC_karyotype.txt"));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        int[][] karyotyLocationMatrix = ves.stream().map(x -> x.stream().mapToInt(Double::intValue).toArray()).toArray(int[][]::new); //converts location data to matrix
        System.out.println("Total populations "+karyotyLocationMatrix.length);
        return karyotyLocationMatrix;
    }

    public static double roundToDecimalPlaces(double value, int decimalPlaces) {
        double scale = Math.pow(10, decimalPlaces);
        return Math.round(value * scale) / scale;
    }



    //Copy of Customized version of MaxDelta() inside  PDEGRID2D.java class
    // Returns both max delta and the grid location
    //run after diffusionADI and before update()

    /**
    public double[] MaxDeltaCustom() {
        double maxDif_cus =deltas[0];
        int maxIndex=0;
        for (int i = 0; i < field.length; i++) {
            if(Math.abs(deltas[i])>maxDif_cus) {
                maxDif_cus = Math.abs(deltas[i]);

                maxIndex=i;
            }
        }
        return new double[] {maxIndex,maxDif_cus};
    }




     public double[] MaxDeltaCustom2(int[][] vssls) {
     List<Integer> excludeList = new ArrayList<>();
     for (int[] row : vssls) {
     int[] maphoodResult = MooreHood(true);
     MapHood(maphoodResult,row[0], row[1]);

     for (int val=0;val<9;val++) {
     excludeList.add(maphoodResult[val]);

     }
     }

     // Return field excluding the indices in excludeList
     return IntStream.range(0, field.length)
     .filter(i -> !excludeList.contains(i))
     .mapToDouble(i -> Math.abs((deltas[i])))
     .toArray();
     }
     */



    /** The next two metods sshould be inserted in PDEGRID2D class. It is used to return the maximum delta
     * at each call of the diffusion model.
    public double[] MaxDeltaCustom2(int[][] vssls) {
        List<Integer> excludeList = new ArrayList<>();

        // Gather exclude indices by applying MapHood
        for (int[] row : vssls) {
            int[] maphoodResult = MooreHood(true);
            MapHood(maphoodResult, row[0], row[1]);

            // Add the result of MapHood to exclude list
            for (int val = 0; val < 9; val++) {
                excludeList.add(maphoodResult[val]);
            }
        }

        // Create a container to hold the max value and its index
        PDEGrid2D.MaxResult maxResult = new PDEGrid2D.MaxResult(Double.NEGATIVE_INFINITY, -1);

        // Stream through the field, excluding the indices in excludeList
        IntStream.range(0, field.length)
                .filter(i -> !excludeList.contains(i))  // Exclude indices
                .forEach(i -> {
                    double absValue = Math.abs(deltas[i]);
                    if (absValue > maxResult.maxValue) {
                        maxResult.maxValue = absValue;
                        maxResult.maxIndex = i;  // Store the original index of the max value
                    }
                });

        // Return both the max value and its original index as a double[]
        return new double[]{maxResult.maxValue, maxResult.maxIndex};
    }

    // Helper class to store the max value and its index
    class MaxResult {
        double maxValue;
        int maxIndex;

        MaxResult(double maxValue, int maxIndex) {
            this.maxValue = maxValue;
            this.maxIndex = maxIndex;
        }
    }

    */

}


