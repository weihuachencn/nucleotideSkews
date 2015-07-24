/*
 * for bacterial codons 
 create on April 29, 2014 --
 */
package preferredcodonanalysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author wchen
 */
public class CodonStats {

    private final HashMap<String, HashMap<Integer, Integer>> hmCodon2DIdx = new HashMap<>();
    private final HashMap<String, String> hmCodon2AA = new HashMap<>();
    private final HashMap<String, ArrayList<String>> hmAA2Codons = new HashMap<>();
    private final ArrayList<String> aas = new ArrayList<>();

    public CodonStats() {
        // codons
        char[] aB1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".toCharArray();
        char[] aB2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".toCharArray();
        char[] aB3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".toCharArray();
        char[] aA  = "FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".toCharArray();

        // get codon to aa hashmap; stop codons are skipped --
        for (int i = 0; i < aB1.length; i++) {
            String aa = aA[i] + "";
            StringBuilder codon = new StringBuilder();
            codon.append(aB1[i]).append(aB2[i]).append(aB3[i]);
            if ( codon.length() == 3) {
                String codonStr = codon.toString();
                hmCodon2AA.put(codonStr, aa);
                
                if( !hmAA2Codons.containsKey(aa) ){
                    hmAA2Codons.put(aa, new ArrayList<String>());
                }
                
                hmAA2Codons.get(aa).add(codonStr);
            }
        }// 
        
        aas.addAll( this.hmAA2Codons.keySet());
        Collections.sort(aas); // 

        // get codon degeneracy index --
        char[] nts = "ATGC".toCharArray();
        for (Map.Entry<String, String> entry : hmCodon2AA.entrySet()) {
            String codon = entry.getKey();
            String aa = entry.getValue();

            /**
             * codon degeneracy index --
             */
            char[] ntsInCodon = codon.toCharArray();

            for (int p = 1; p <= 3; p++) {
                int dindex = 0;
                for (char nt : nts) {
                    StringBuilder newcodon = new StringBuilder();
                    if (p == 1) {
                        newcodon.append(nt).append(ntsInCodon[1]).append(ntsInCodon[2]);
                    } else if (p == 2) {
                        newcodon.append(ntsInCodon[0]).append(nt).append(ntsInCodon[2]);
                    } else if (p == 3) {
                        newcodon.append(ntsInCodon[0]).append(ntsInCodon[1]).append(nt);
                    }

                    String ncodon = newcodon.toString();
                    String newaa = hmCodon2AA.containsKey(ncodon) ? hmCodon2AA.get(ncodon) : "X";
                    if (newaa.equalsIgnoreCase(aa)) {
                        dindex++;
                    }
                } // end of for each nucleotide

                // -- --
                if (!hmCodon2DIdx.containsKey(codon)) {
                    hmCodon2DIdx.put(codon, new HashMap<Integer, Integer>());
                }
                hmCodon2DIdx.get(codon).put(p, dindex);
            }// iterate each codon position
        }// iterate each codon
    } // the constructor 
    
    public HashMap<String, HashMap<Integer, Integer>>  getCodon2DIdx(){
        return this.hmCodon2DIdx;
    }
    
    public HashMap<Integer,Integer> getDIdx( String codon ){
        return this.hmCodon2DIdx.get(codon);
    }
    
    public HashMap<String, String> getCodon2AA(){
        return this.hmCodon2AA;
    }
    
    public HashMap<String, ArrayList<String>> getAA2Codons(){
        return this.hmAA2Codons;
    }
    
    public ArrayList<String> getCodons4AAminoAcid(String aa){
        return this.hmAA2Codons.get(aa);
    }
    
    public ArrayList<String> getAAs(){
        return this.aas;
    }
    
    public String getAAbyCodon( String codon ){
        return this.hmCodon2AA.get(codon);
    }
    
}
