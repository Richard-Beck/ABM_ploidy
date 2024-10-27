package onLatticeCA_jar;
import java.util.*;

public class homeostasisTests {
    // We may wish to pass in the following 2 parameters via params file.
    int winsize =20;
    double criticalValue = 1.;

    boolean homeostasisReached = false;
    boolean homeostasisAchieved=false;
    List<Double> o2List;
    List<Double> nCellsList;

    homeostasisTests(){
        o2List = new ArrayList<>();
        nCellsList = new ArrayList<>();
    }
    public boolean test(double o2, int ncells){
        if(homeostasisReached) return true;
        o2List.add(0,o2);
        nCellsList.add(0,(double)ncells);
        if(o2List.size()<2*winsize) return false;
        if(o2List.size()>2*winsize) o2List.remove(o2List.size()-1);
        if(nCellsList.size()>2*winsize) nCellsList.remove(nCellsList.size()-1);
         //homeostasisReached = ks_test(o2List)&ks_test(nCellsList);
         homeostasisReached=checkHomeostasis(o2List,winsize,1e-4)&checkHomeostasis(nCellsList,winsize,100);


        if(homeostasisReached){
          //  System.out.println("homeostasis has been reached!");
        }
        return homeostasisReached;
    }

    private boolean ks_test(List<Double> list){
        double[] sample1 = new double[winsize];
        double[] sample2 = new double[winsize];
        for(int i = 0; i<winsize;i++) sample1[i] = list.get(i);
        for(int i = 0; i<winsize;i++) sample2[i] = list.get(i+winsize);
        List<Double> allValues = getSortedUniqueValues(list);
        double maxDiff = 0.;
        for (Double value : allValues) {
            double delta = ecdf(sample1, value) - ecdf(sample2, value);
            maxDiff = Math.max(Math.abs(delta), maxDiff);
        }
        double testStat = criticalValue*Math.sqrt((double) (sample1.length + sample2.length) /(sample1.length*sample2.length));
        return maxDiff<testStat;
    }

    private List<Double> getSortedUniqueValues(List<Double> list){
        // hashset auto filters unique values
        Set<Double> uniqueSet = new HashSet<>(list);
        List<Double> uniqueList = new ArrayList<>(uniqueSet);
        Collections.sort(uniqueList);
        return(uniqueList);
    }

    private double ecdf(double[] data, double value){
        double sum=0.;
        for (double datum : data) {
            if (datum <= value) sum += 1.;
        }
        return sum/data.length;
    }

/** Sliding window test, check the overall difference between the average of two consecutive windows in O2 values and compares
 * it to a set tolerance  */


    public static boolean checkHomeostasis(List<Double> o_data, int windowSize, double tolerance) {
        // Convert List to an array for easier access
        double[] data = new double[o_data.size()];
        for (int i = 0; i < o_data.size(); i++) {
            data[i] = o_data.get(i);
        }

        // Make sure the data length is at least 2 times the window size
        if (data.length < windowSize * 2) {
            throw new IllegalArgumentException("Data array is too small for the given window size.");
        }

        // Get the last two windows
        int secondWindowStart = data.length - windowSize; // Start of the second window
        int firstWindowStart = secondWindowStart - windowSize; // Start of the first window

        // Calculate the average value in the first window [firstWindowStart, secondWindowStart)
        double avgWindow1 = calculateAverage(data, firstWindowStart, secondWindowStart);

        // Calculate the average value in the second window [secondWindowStart, data.length)
        double avgWindow2 = calculateAverage(data, secondWindowStart, data.length);

        // Compute the absolute change between the two consecutive windows
        double avgChange = Math.abs(avgWindow2 - avgWindow1);

        // If the change is greater than the tolerance, homeostasis has not been reached
        return avgChange <= tolerance;
    }



    private static double calculateAverage(double[] data, int start, int end) {
        double sum = 0;
        for (int i = start; i < end; i++) {
            sum += data[i];
        }
        return sum / (end - start);
    }

}
