// ========================================================================
// Definition of a single cell.
// ========================================================================

package onLatticeCA_jar;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.Gui.UIGrid;
import HAL.Util;
import HAL.Rand;

import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;


import HAL.GridsAndAgents.PDEGrid2D;

// ========================================================================
// Class to model individual cells

/**
 * This is the main AGENT class;
 * Manages cell division and assignment of various agent attributes, eg. division rate, death rate etc
 * Contain a draw method which visualizes agent in a UIGrid
 */

public class  Cell extends AgentSQ2Dunstackable<OnLatticeGrid> {
    double divisionRate;


    double movementRate;
    double deathRate;
    int resistance; // 0 = sensitive, 1 = resistant
    public double missegregationRate;//=.125; // Rate at which cells missegregate
    public boolean cancerous;

    double Vo; // uptake mmol/L/S
    double ko; // michaelis const. mmol/L

    final public int[] karyotype = new int[22];

    public int ChromToChange;

    private int copyNumberChange;

    public int hash;
    helperMethods helper = new helperMethods();
    public double[] consumption = new double[3];

    public int[][] bloodVessels;
    // -----------------------------------
    public void init(int res, int[] kary){
        resistance = res;
        for(int i = 0; i<karyotype.length; i++){
            karyotype[i] = kary[i];
        }
        ko = Params.consumption_MM_constant;
        hash = Arrays.hashCode(karyotype);   // unique identity generated for each cell based on karyotype
        deathRate = Params.natural_death_rate;
        if (Util.ArrayMean(karyotype) == 2) {
            cancerous = false;
            divisionRate = Params.natural_growth_rate;
            Vo = Params.normal_oxygen_consumption_rate;
        } else {
            cancerous = true;
            divisionRate = Params.cancer_cell_growth_rate;
            Vo = Params.cancer_oxygen_consumption_rate;
        }

        setMissegregationRate(karyotype, cancerous);

        // I'm not sure if this resistance actually does anything, but better leave it in
        // for now
        if (resistance == 0) {
            movementRate = G.movementRate_S;
        } else {
            movementRate = G.movementRate_R;
        }
    }
    public void init(Cell parent){
        // copy over parameters from parent. Would be better if we could just copy
        // the parent then make any modifications after if necessary. We would have to check
        // the newAgent constructor to see if it is doing any other tasks that we would
        // need to assume responsibility for / RJB
        cancerous = parent.cancerous; Vo=parent.Vo; ko=parent.ko;
        movementRate = parent.movementRate; resistance = parent.resistance;
        missegregationRate = parent.missegregationRate;
        //copy karyotype from parent
        for(int i = 0; i<karyotype.length; i++){
            karyotype[i] = parent.karyotype[i];
        }
    }
    boolean HasSpace() {
        if(cancerous){
            return MapEmptyHood(G.hood) > 0;
        }else{
            return MapEmptyHood(G.hood) >1;
        }
    }

    private void missegregate(Cell daughter){

        cancerous = true;
        daughter.cancerous=true;

        ChromToChange = G.rn.Int(karyotype.length);
        copyNumberChange =G.rn.Int(karyotype[ChromToChange]);

        karyotype[ChromToChange] = 2 * karyotype[ChromToChange] - copyNumberChange;
        daughter.karyotype[ChromToChange] = copyNumberChange;

        daughter.setMissegregationRate(daughter.karyotype,daughter.cancerous);//
        setMissegregationRate(karyotype, cancerous);
        G.missegregationCounter += 1;
    }

    void make_cancerous(){
        cancerous= true;
        Vo = Params.cancer_oxygen_consumption_rate;
        divisionRate = Params.cancer_cell_growth_rate;
    }
    void Divide () {
        // How does this code behave when the neighborhood is full?
        // OR is it the case that we can never get here when the neighborhood is full...
        int nOpts=MapEmptyHood(G.hood);
        int d_opts= G.rn.Int(nOpts);
        int iDaughter = G.hood[d_opts];

        // abort if chosen daughter location overlaps with vessel
        if(helper.IsMember(iDaughter,G.resources.vessel_indices)) return;

        Cell daughter = G.NewAgentSQ(iDaughter);
        // copy parent parameters to daughter
        daughter.init(this);
        if (G.rn.Double() < missegregationRate & G.homeostasisReached){
            missegregate(daughter);
            make_cancerous();
            daughter.make_cancerous();

        } else G.RegularDivisionCounter += 1;

        hash = Arrays.hashCode(karyotype);
        daughter.hash = Arrays.hashCode(daughter.karyotype);

        //ISSUES:
        // This 5 should not be hardcoded here. Better places include but are not limited
        // to the parameter file, or as a static int at the top of a class definition.
        // Isn't it possible the parent is also above the threshold? Why isn't that checked??
        // Seems there is another check inside the stepcells function, which should catch the parent.
        // Probably should delete the check below.
        if (helper.isAboveThreshold(5, daughter.karyotype)) {
            daughter.Dispose();
            G.DeadCellCounter +=1;
            // what does this do.
            G.cellCountsArr[daughter.resistance] -= 1;
        } else {
            // seems a little odd to draw the daughter now.
            daughter.Draw();
        }

    }
    // -----------------------------------
    boolean Move(){
        boolean successfulMoveB = false;
        int nOpts=MapEmptyHood(G.hood);// identify the empty spots in the cell's neighbourhood (if there are any).
        if(nOpts>0){
           int Opts_int=G.rn.Int(nOpts);//G.rn.Int(nOpts);
           G.mainRnNumbers.add((double)Opts_int);
            int iDestination = G.hood[Opts_int];

            MoveSafeSQ(G.ItoX(iDestination), G.ItoY(iDestination));
            successfulMoveB = true;
        }
        return successfulMoveB;
    }
  public void setKaryotype(int[] Array){
        for(int i=0;i<Array.length;i++){
            karyotype[i]=Array[i];
        }

    }
    public void setMissegregationRate(int [] karyotype, boolean type){
        //Use array sum to inform missegregation rate
        if(cancerous){
            missegregationRate =(Util.ArraySum(karyotype)*Params.cancer_cell_mis_seg_rate);
        }else{
            missegregationRate =Params.normal_cell_mis_seg_rate;//(Util.ArraySum(karyotype)*Params.natural_mis_seg_rate)/(Params.misRateWeight);

        }
    }


    public void setDivisionRate(){ //concentration dependent division rate
        double O2concentration=G.resources.pdegrid2d.Get(Isq());
        if(cancerous){
            divisionRate=(O2concentration * Params.cancer_cell_growth_rate) /
                    (Util.ArrayMean(karyotype) * Params.growth_michaelis_menten_constant + O2concentration);
        }else {
            divisionRate = (O2concentration * Params.natural_growth_rate) /
                    ((Util.ArrayMean(karyotype) * Params.growth_michaelis_menten_constant + O2concentration));
        }
    }

    public void setDeathRate(){ //concentration dependent death rate
        double O2concentration=G.resources.pdegrid2d.Get(Isq());
        double death=(Params.natural_death_rate)*(1-O2concentration/ ( O2concentration + Params.death_michaelis_menten_constant));
        if(cancerous) {
            deathRate = death/20;
        }else{
            deathRate=death;
        }
    }


    public  void Consumption( double timescalar){
        double O= G.resources.pdegrid2d.Get(Isq());
        // oxygen consumption follows Michaelis Menten kinetics:
        double consumption= -Vo*O/(O+ko);
        G.resources.pdegrid2d.Add(Isq(), consumption*timescalar);
    }

    public void DivideOrDie(){
        setDeathRate();
        setDivisionRate();
        double r = G.rn.Double();
        double r2 = G.rn.Double();
        if(r<r2){
            if (r < divisionRate && HasSpace()) {// else if (r< totPropensity) {
                Divide();
                G.divisionCounter++;
            }
        }else{
            if(r2<deathRate){
                Dispose();
                G.vis.SetPix(Isq(),G.BLACK);
                G.DeadCellCounter ++;
            }
        }
    }







    // -----------------------------------
    //Draw cells in different clors based on karyotype
    // For colours from graphs: Util.RGB256(17, 91, 28): Util.RGB256(132, 18, 10)
    // For red/blue combination: Util.RGB256(26, 133, 255): Util.RGB256(212, 17, 89))
    void Draw(){
       // @TODO MODIFY CODE TO DRAW EACH KARYOTYPE
        //G.vis.SetPix(Isq(), Util.HeatMapJet(hash,-999999999,999999999));//sets a single pixel
        G.vis.SetPix(Isq(), (!cancerous) ? Util.RGB256(0, 0, 255) : Util.RGB256(216, 27, 96));
//        G.vis.SetPix(Isq(), (resistance==0)? Util.CategorialColor(13): Util.CategorialColor(11));//sets a single pixel

       // G.vis.SetPix(Isq(), Util.CategorialColor(karyotype));//sets a single pixel
//        G.vis.SetPix(Isq(), (resistance==0)? Util.CategorialColor(13): Util.CategorialColor(11));//sets a single pixel


        //G.vis.SetPix(Isq(), (resistance==0)? Util.RGB256(117, 197, 114): Util.RGB256(216, 27, 96));//sets a single pixel
//        G.vis.SetPix(Isq(), (resistance==0)? Util.CategorialColor(13): Util.CategorialColor(11));//sets a single pixel
    }
   //void DrawPloidy(){
    //G.vis.SetPix(Isq(), (==0)? Util.RGB256(51, 153, 255): Util.RGB256(255, 204, 0));//sets a single pixel
//        G.vis.SetPix(Isq(), (resistance==0)? Util.CategorialColor(13): Util.CategorialColor(11));//sets a single pixel
    //}

}


