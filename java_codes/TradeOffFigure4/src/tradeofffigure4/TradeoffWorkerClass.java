/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tradeofffigure4;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import org.apache.commons.lang3.StringUtils;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.transcription.TranscriptionEngine;

/**
 * see the main class for last modified dates and information ...
 * @author wchen
 */
public class TradeoffWorkerClass implements Callable<String> {

    /**
     * global variables
     */
    // transcription enging / codon table 
    private final TranscriptionEngine engine = new TranscriptionEngine.Builder().table(11).translateNCodons(false).build();
    // cost of aas --

    private float gc;
    private int bias, tasks;

    private int valid_codons = 1000000;

    private HashMap<String, HashMap<Integer, Integer>> hmCodonStates;
    private ArrayList<HashMap<String, Float>> arAA2Cost;
    private String skews = "all"; // Sep 21, 2014; all = calculate both at and gc skews; atonly, gconly = calculate only at or gc skews

    public TradeoffWorkerClass(float gc, int bias, HashMap<String, HashMap<Integer, Integer>> hmCodonStates, 
            ArrayList<HashMap<String, Float>> hmAA2Cost, int tasks, String skews) {
        this.bias = bias;
        this.gc = gc;
        this.arAA2Cost = hmAA2Cost;
        this.hmCodonStates = hmCodonStates;
        this.tasks = tasks;
    }

    // Sep 21, 2014; added options atonly or gc only ...
    @Override
    public String call() throws Exception {
        int gs = !skews.equalsIgnoreCase( "atonly" ) ? 
                (int) (gc / 100 * valid_codons * bias / 100) : 
                (int) (gc / 100 * valid_codons * 50 / 100);;
        int cs = !skews.equalsIgnoreCase( "atonly" ) ? 
                (int) (gc / 100 * valid_codons * (100 - bias) / 100) : 
                (int) (gc / 100 * valid_codons * ( 100 - 50) / 100);
        int as = !skews.equalsIgnoreCase( "gconly" ) ? 
                (int) ((100 - gc) / 100 * valid_codons * bias / 100) : 
                (int) (( 100 - gc ) / 100 * valid_codons * 50 / 100);
        int ts = !skews.equalsIgnoreCase( "gconly" ) ? 
                (int) ((100 - gc) / 100 * valid_codons * (100 - bias) / 100) : 
                (int) (( 100 - gc ) / 100 * valid_codons * ( 100 - 50) / 100);

        // generate ranom coding sequence
        StringBuilder sb = new StringBuilder();
        sb.append(StringUtils.repeat("A", as * 2)).
                append(StringUtils.repeat("T", ts * 2)).
                append(StringUtils.repeat("G", gs * 2)).
                append(StringUtils.repeat("C", cs * 2));

        ArrayList<Character> al = new ArrayList<>();
        for (char c : sb.toString().toCharArray()) {
            al.add(c);
        }
        Collections.shuffle(al); // shuffle!!
        String seq = StringUtils.join(al, "");

        DNASequence seqobj = new DNASequence(seq, DNACompoundSet.getDNACompoundSet());
        ProteinSequence proSeq = seqobj.getRNASequence().getProteinSequence(engine);
        String proSeqStr = proSeq.getSequenceAsString();

        /**
         * calculate: skews, cost of AAs
         */
        HashMap<String, Integer> hmAA2Counts = new HashMap<>();
        for (AminoAcidCompound aac : AminoAcidCompoundSet.getAminoAcidCompoundSet().getAllCompounds()) {
            String aaShortName = aac.getShortName();
            hmAA2Counts.put(aaShortName, 0); // update number
        }// foreach AA

        int current_valid_codon = 0;
        HashMap<String, Integer> hmCodon2Counts = new HashMap<>();
        for (String codon : hmCodonStates.keySet()) {
            hmCodon2Counts.put(codon, 0);
        }

        for (int idx = 0; idx < proSeqStr.length(); idx++) {
            String aa = proSeqStr.charAt(idx) + "";
            if (!aa.equalsIgnoreCase("*") && arAA2Cost.get(0).containsKey(aa)) {
                current_valid_codon++;

                // count AAs 
                hmAA2Counts.put(aa, hmAA2Counts.get(aa) + 1);

                // count codon
                String codon = seq.substring(idx * 3, idx * 3 + 3);
                hmCodon2Counts.put(codon, hmCodon2Counts.get(codon) + 1);
            }// if valid 

            if (current_valid_codon >= this.valid_codons) {
                break;
            }
        }//

        ArrayList<Float> skews = skewsCacl(hmCodon2Counts);
        ArrayList<Float> meanAACost = meanAACostCacl(hmAA2Counts);

        /**
         * final results
         */
        StringBuilder sb2 = new StringBuilder();
        sb2.append(gc);
        sb2.append("\t").append(bias);
        sb2.append("\t").append((bias - (100 - bias)) / 100f); // skews 

        for (float f : skews) {
            sb2.append("\t").append(f);
        }
        
        for( float f : meanAACost ){
            sb2.append("\t").append(f);
        }
        sb2.append("\n"); // also append a new line 

        if (tasks % 100 == 0) {
            System.err.println("----------------------------------------------------------\n\t"
                    + tasks + " tasks done \n"
                    + "----------------------------------------------------------");
        }
        return sb2.toString();
    }

    /**
     *
     * @param hmCodon2Counts
     * @return
     */
    private ArrayList<Float> skewsCacl(HashMap<String, Integer> hmCodon2Counts) {
        HashMap<Integer, HashMap<String, Integer>> hmIdx2Nt2Count = new HashMap<>();

        for (Map.Entry<String, Integer> entry : hmCodon2Counts.entrySet()) {
            String codon = entry.getKey();
            int count = entry.getValue();
            if (hmCodonStates.containsKey(codon)) {
                for (int i = 1; i <= 3; i++) {
                    String letter = codon.charAt(i - 1) + "";
                    int didx = hmCodonStates.get(codon).get(i);
                    if (!hmIdx2Nt2Count.containsKey(didx)) {
                        hmIdx2Nt2Count.put(didx, new HashMap<String, Integer>());
                    }
                    if (!hmIdx2Nt2Count.get(didx).containsKey(letter)) {
                        hmIdx2Nt2Count.get(didx).put(letter, 0);
                    }
                    hmIdx2Nt2Count.get(didx).put(letter, hmIdx2Nt2Count.get(didx).get(letter) + count);
                }// for each codon position
            }// if codon is valid
        }// foreach codon

        ArrayList<Float> first = calcSkew(hmIdx2Nt2Count, 1);
        ArrayList<Float> second = calcSkew(hmIdx2Nt2Count, 2);
        ArrayList<Float> fourth = calcSkew(hmIdx2Nt2Count, 4);

        ArrayList<Float> skews = new ArrayList<>();
        skews.addAll(first);
        skews.addAll(second);
        skews.addAll(fourth);

        return skews;
    }// end of skewsCacl

    /**
     * -- a private function --
     * @param hmAA2Counts; $hash{ $aa } = $number_of_aas ;
     * @return 
     */
    private ArrayList<Float> meanAACostCacl(HashMap<String, Integer> hmAA2Counts) {

        ArrayList<Float> meanAAs = new ArrayList<>();

        for (HashMap<String, Float> aa2cost : arAA2Cost) {
            float total_enery = 0f;
            int total_aas = 0;
            for (Map.Entry<String, Integer> entry : hmAA2Counts.entrySet()) {
                String aa = entry.getKey();
                int count = entry.getValue();
                if (aa2cost.containsKey(aa) && count > 0) {
                    total_aas += count;
                    total_enery += count * aa2cost.get(aa);
                }
            }
            
            float mean = total_aas == 0 ? 0 : total_enery / total_aas;
            
            //
            meanAAs.add(mean);
        }
        return meanAAs;
    }// end of meanAACostCacl

    private ArrayList<Float> calcSkew(HashMap<Integer, HashMap<String, Integer>> hmIdx2Nt2Count, int i) {
        ArrayList<Float> skews = new ArrayList<>();
        if (hmIdx2Nt2Count.containsKey(i)) {
            HashMap<String, Integer> nt2count = hmIdx2Nt2Count.get(i);
            int a = nt2count.containsKey("A") ? nt2count.get("A") : 0;
            int t = nt2count.containsKey("T") ? nt2count.get("T") : 0;
            int g = nt2count.containsKey("G") ? nt2count.get("G") : 0;
            int c = nt2count.containsKey("C") ? nt2count.get("C") : 0;
            int at = a + t;
            int gc = g + c;

            // calculate skews 
            float atskew = at == 0 ? 0f : (a - t) / (float) at;
            float gcskew = gc == 0 ? 0f : (g - c) / (float) gc;
            skews.add(atskew);
            skews.add(gcskew);
        } else {
            skews = new ArrayList<>(Arrays.asList(0f, 0f));
        }
        return skews;
    }

}
