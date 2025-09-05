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
    void missegregate(int N, int X, Cell daughter){
        final int nChr = karyotype.length;

        // snapshot of parent's pre-division counts to define event probabilities
        int[] base = java.util.Arrays.copyOf(karyotype, nChr);

        // cumulative sums over base to map [0..N-1] -> chromosome index
        int[] cum = new int[nChr + 1];
        cum[0] = 0;
        for (int j = 0; j < nChr; j++) cum[j + 1] = cum[j] + base[j];

        // X missegregation events
        for (int e = 0; e < X; e++) {
            int pick = G.rn.Int(N);            // 0..N-1
            // map to chromosome j (small nChr => linear scan is fine)
            int j = 0; while (pick >= cum[j + 1]) j++;

            // choose direction randomly (which daughter gains)
            boolean daughterGains = G.rn.Bool();

            if (daughterGains) {
                // parent -> daughter if parent has something to give
                if (karyotype[j] > 0) {
                    karyotype[j]      -= 1;
                    daughter.karyotype[j] += 1;
                } else if (daughter.karyotype[j] > 0) {
                    // flip direction to avoid negative
                    daughter.karyotype[j] -= 1;
                    karyotype[j]      += 1;
                }
            } else {
                // daughter -> parent if daughter has something to give
                if (daughter.karyotype[j] > 0) {
                    daughter.karyotype[j] -= 1;
                    karyotype[j]      += 1;
                } else if (karyotype[j] > 0) {
                    // flip direction
                    karyotype[j]      -= 1;
                    daughter.karyotype[j] += 1;
                }
            }
        }

        daughter.setMissegregationRate(daughter.karyotype, daughter.cancerous);
        setMissegregationRate(karyotype, cancerous);
        G.missegregationCounter += 1;
    }   


    void make_cancerous(){
        cancerous= true;
        Vo = Params.cancer_oxygen_consumption_rate;
        divisionRate = Params.cancer_cell_growth_rate;
    }
    void Divide() {
        int nOpts = MapEmptyHood(G.hood);
        if (nOpts == 0) return;

        int iDaughter = G.hood[G.rn.Int(nOpts)];
        if (helper.IsMember(iDaughter, G.resources.vessel_indices)) return;

        Cell daughter = G.NewAgentSQ(iDaughter);
        daughter.init(this);

        // seed exactly one cancer at/after homeostasis
        if (G.homeostasisReached && !G.cancerExists) {
            make_cancerous();
            G.cancerExists = true;
        }
        if (this.cancerous) daughter.make_cancerous();

        if (this.cancerous && G.cancerExists) {
            // 1) sum copies
            int N = 0; for (int v : karyotype) N += v;
            if (N > 0) {
                // 2) one binomial draw using RAND.java
                int X = G.rn.Binomial(N, missegregationRate);
                // 3) allocate X events among N copies
                if (X > 0) missegregate(N, X, daughter);
            }
        } else {
            G.RegularDivisionCounter += 1;
        }

        // refresh hashes 
        hash = Arrays.hashCode(karyotype);
        daughter.hash = Arrays.hashCode(daughter.karyotype);

        this.Draw();
        daughter.Draw();
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
            missegregationRate =Params.cancer_cell_mis_seg_rate;
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
        if (!cancerous) {
            G.vis.SetPix(Isq(), HAL.Util.RGB256(0, 0, 255));   // normal: blue
        } else {
            // cancer: color by karyotype identity
            G.vis.SetPix(Isq(), G.colorForKaryotypeHash(hash));
        }
    }

}


