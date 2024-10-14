package onLatticeCA_jar;
import java.util.*;

public class homeostasisTests {
    // We may wish to pass in the following 2 parameters via params file.
    int winsize =20;
    double criticalValue = 1.;

    boolean homeostasisReached = false;
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
        homeostasisReached = ks_test(o2List)&ks_test(nCellsList);
        if(homeostasisReached){
            System.out.println("homeostasis has been reached!");
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
}
