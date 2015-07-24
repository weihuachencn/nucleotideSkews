/*
 * this is a single thread version...
 * 
 */
package preferredcodonsh;

import com.google.common.base.Joiner;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

/**
 *
 * @author wchen
 * last modified : Nov 4, 2014 
 * this is the complete version!!
 */
public class PreferredCodonSH {

    private final CodonStats cs = new CodonStats(); // codon stats --
    private final ArrayList<String> allCodons = new ArrayList<String>();
    private HashMap<String, String> codon2aa;
    
    private final HashMap<String, HashMap<String, Integer>> aa2gene2count = new HashMap<>();

    public PreferredCodonSH() {
        String codonsString = "TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG";
        for (String codon : codonsString.split("\\s+")) {
            allCodons.add(codon);
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParseException, IOException, Exception {
        new PreferredCodonSH().run(args);
    }

    private void run(String[] args) throws ParseException, IOException, InterruptedException, Exception {
        /**
         * command line parser using CLI.commons
         */
        // create the command line parser
        CommandLineParser parser = new BasicParser();

        // create the Options
        Options o = new Options();

        o.addOption("i", "codingseq-list-file", true, "accession number to coding sequence list; mandatory");
        o.addOption("s", "seqcount-exe", true, "path to the executable of SeqCount; mandatory");
        o.addOption("e", "encprime-exe", true, "path to the executable of ENCprime; mandatory");
        o.addOption("o", "output-file", true, "output file; mandatory");
        o.addOption("d", "debug", false, "debug mode, only do first 5 jobs; optional");
        
        /**
         * parse options and check if all required options are there
         */
        CommandLine cmd = parser.parse(o, args);

        if (!cmd.hasOption("c") || !cmd.hasOption("o") ||  !cmd.hasOption("s") || !cmd.hasOption("e")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("this.jar", o);
            System.exit(1);
        }

        if (cmd.hasOption("d")) {
            System.err.println("\n\tdebug mode, only the first 5 jobs will be submitted\n\n");
        }
        
        /**
         * read command line parameters
         */
        String acc2seqfile = cmd.getOptionValue("i");
        String outfile = cmd.getOptionValue("o");

        // acc to seq files 
        HashMap<String, String> acc2CDSfile = readListFile( acc2seqfile );

        // get exe files --
        String SeqCountExe = cmd.getOptionValue("s");
        String ENCprimeExe = cmd.getOptionValue("e");
        
        /**
         * prepare output file
         */
        BufferedWriter out;
        out = new BufferedWriter(new FileWriter(outfile));
        String joined = Joiner.on("\t").join( this.cs.getAAs() );
        out.write(joined);
        out.newLine();
        
        /**
         * 
         */
        int tasks = 0;


        for (Map.Entry<String, String> entry : acc2CDSfile.entrySet()) {
            String acc = entry.getKey();
            String fastafile = entry.getValue();
            tasks++;

            if (tasks % 100 == 0) {
                System.err.println("----------------------------------------------------------\n\t"
                        + tasks + " tasks done \n"
                        + "----------------------------------------------------------");
            }

            this.codon2aa = cs.getCodon2AA();

            // step 1. link file to current dir --
            Process p1 = Runtime.getRuntime().exec("ln -s " + fastafile + " .");
            p1.waitFor(); // wait for the process to finish ...

            // step 2. count how many genes in the fasta file
            ArrayList<String> geneNames = readFastaForGeneNames( fastafile );
            int nseqs = geneNames.size();

            System.err.println("==> task: " + tasks + " ------- " + acc + " " + nseqs + " genes ----------");
            
            String codCnt = acc + ".codcnt";
            String atgcCnt = acc + ".acgtcnt";
            String resultFile = acc + ".results";
            
            // step 3. run SeqCount and ENCprime
            Process p3a = Runtime.getRuntime().exec(SeqCountExe + " -c " + acc + " " + nseqs);
            p3a.waitFor(); // wait for the process to finish ...

            Process p3b = Runtime.getRuntime().exec(SeqCountExe + " -n " + acc + " " + nseqs);
            p3b.waitFor(); // wait for the process to finish ...

            String exe = ENCprimeExe + " " + codCnt + " " + atgcCnt + " 11 " + resultFile + " 0 -q";
            Process p3c = Runtime.getRuntime().exec(exe);
            p3c.waitFor(); // wait for the process to finish ...

            // step 4, load data and do analyses ; 
            // but first, check if all file exists ... 
            File conCntFile = new File(codCnt);
            File resultFileFile = new File(resultFile);
            if (conCntFile.exists() && resultFileFile.exists()) {
                /**
                 * load resultFile --
                 */
                HashMap<String, ArrayList<Double>> gene2nc = loadResultFile(resultFile);
                System.err.println("\t --> load gene2nc : " + gene2nc.size());

                /**
                 * load codon count file -- $hash{ $codon }{ $gene } = $number
                 *
                 * $hash{ $aa }{ $gene } = $total_count ...
                 */
                HashMap<String, HashMap<String, Integer>> hmCodon2Gene2Count = loadCodonCountFile(conCntFile);
                System.err.println("\t --> load codonCountFile : " + hmCodon2Gene2Count.size());

                /**
                 * clean up...
                 */
                conCntFile.delete();
                new File(acc + ".codfreq").delete();
                new File(acc + ".acgtcnt").delete();
                new File(acc + ".acgtfreq").delete();
                new File(acc).delete();
                resultFileFile.delete();

                // analyses 
                StringBuilder sb = new StringBuilder();
                sb.append(acc);

                for (String aa : this.cs.getAAs()) {
                    String preferredCodon = "none";
                    double current_rho = 10;

                    if (this.aa2gene2count.containsKey(aa)) {
                        ArrayList<String> codons = cs.getCodons4AAminoAcid(aa);
                        int n = codons.size();
                    //System.err.println( "\t\t" + aa + " - > " + n );

                        HashMap<String, Integer> gene2AAFamCount = this.aa2gene2count.get(aa);

                        for (String codon : codons) {
                            ArrayList<Double> nc = new ArrayList<Double>(), ncp = new ArrayList<>(), y = new ArrayList<>();
                            
                            for (String gene : geneNames) {
                                int codoncount = hmCodon2Gene2Count.containsKey(codon) && hmCodon2Gene2Count.get(codon).containsKey(gene)
                                        ? hmCodon2Gene2Count.get(codon).get(gene) : 0;

                                int codonFamCount = gene2AAFamCount.containsKey(gene) ? gene2AAFamCount.get(gene) : 0;
                                //System.err.println( "\t\t-->" + codoncount + "; fam count: " + codonFamCount );

                                if (gene2nc.containsKey(gene) && codonFamCount >= 10) {
                                    ArrayList<Double> ncs = gene2nc.get(gene);
                                    nc.add(ncs.get(0));
                                    ncp.add(ncs.get(1));
                                    y.add((double) codoncount); // conver it to double 
                                }
                            }// for each gene

                            /**
                             * correlate nc to y, and ncp to y
                             */
                            if (nc.size() >= 5) {
                                /**
                                 * correlate nc with codon frequency / count
                                 */
                                RealMatrix dnc = makeRealMatrix(nc, y);
                                SpearmansCorrelation spear1 = new SpearmansCorrelation(dnc);
                                PearsonsCorrelation pearson1 = spear1.getRankCorrelation();
                                double rho1 = spear1.getCorrelationMatrix().getEntry(0, 1);
                                double pvalue1 = pearson1.getCorrelationPValues().getEntry(0, 1);

                                if (pvalue1 < 0.05 / n && rho1 < current_rho) {
                                    current_rho = rho1;
                                    preferredCodon = codon;
                                }

                                /**
                                 * correlate ncp with codon frequency / count
                                 */
                                RealMatrix dncp = makeRealMatrix(ncp, y);
                                SpearmansCorrelation spear2 = new SpearmansCorrelation(dncp);
                                PearsonsCorrelation pearson2 = spear2.getRankCorrelation();
                                double rho2 = spear2.getCorrelationMatrix().getEntry(0, 1);
                                double pvalue2 = pearson2.getCorrelationPValues().getEntry(0, 1);

                                if (pvalue2 < 0.05 / n && rho2 < current_rho) {
                                    current_rho = rho2;
                                    preferredCodon = codon;
                                }
                            }
                        } // for each codon
                    }// if exists codon fam 

                    sb.append("\t").append(preferredCodon);
                }// for each aa 
                sb.append("\n");
                out.write(sb.toString()); // the results contain new lines already ; Dec 3, 2013 --
                System.err.println();
            } else {
                System.err.println("something wrong when processing : " + acc);
            }

            if (cmd.hasOption("d") && tasks >= 5) {
                System.err.println("\n\tstop submitting jobs because you're running in debug mode!!\n\n");
                break; // break the loop !!
            }
        } // for each; to submit jobs --
        out.close();
    }

    /**
     * @param file contains a list of accession to value
     * @return a hashmap $hash{ $acc } = $location_of_file
     * @throws FileNotFoundException
     * @throws IOException
     */
    public HashMap<String, String> readListFile(String file) throws FileNotFoundException, IOException {
        HashMap<String, String> acc2file = new HashMap<String, String>();
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            String[] parts = line.trim().split("\t");
            if (parts.length >= 2) {
                acc2file.put(parts[0].trim(), parts[1].trim()); // also trim all the strings 
            }
        }
        br.close();
        return acc2file;
    }

    /**
     * output: $hash{ $gene } = [$nc, $ncp]; double
     */
    private HashMap<String, ArrayList<Double>> loadResultFile(String infile) throws FileNotFoundException, IOException {
        HashMap<String, ArrayList<Double>> hm = new HashMap<>();

        int lines = 0;
        // file1 --
        BufferedReader br = new BufferedReader(new FileReader(infile));
        String line;
        while ((line = br.readLine()) != null) {
            lines++;

            if (lines > 1) { // from the second line on ... 
                String[] parts = line.trim().split(":*\\s+");
                if (parts.length >= 3) {
                    String gene = parts[0].trim();
                    double nc = Double.parseDouble(parts[1].trim());
                    double ncp = Double.parseDouble(parts[2].trim());

                    ArrayList<Double> ns = new ArrayList<>();
                    ns.add(nc);
                    ns.add(ncp);

                    hm.put(gene, ns);
                }
            } // skip title line --
        }
        br.close();
        return hm;
    } // 

    private HashMap<String, HashMap<String, Integer>> loadCodonCountFile(File conCntFile) throws FileNotFoundException, IOException {
        HashMap<String, HashMap<String, Integer>> hm = new HashMap<>();

        int lines = 0;
        // file1 --
        BufferedReader br = new BufferedReader(new FileReader(conCntFile));
        String line;
        while ((line = br.readLine()) != null) {
            lines++;

            // from the fourth line on ....
            if (lines >= 4) {
                String[] parts = line.trim().split(">");
                if (parts.length >= 2) {
                    String gene = parts[0].trim();
                    String[] counts = parts[1].trim().split("\\s+");

                    if (counts.length >= 64 && !gene.equalsIgnoreCase("Totals")) {
                        for (int idx = 0; idx < counts.length; idx++) {
                            String codon = this.allCodons.get(idx);
                            int count = Integer.parseInt(counts[idx]);
                            if (!hm.containsKey(codon)) {
                                hm.put(codon, new HashMap<String, Integer>());
                            }
                            hm.get(codon).put(gene, count);

                            /**
                             * aa to gene to count
                             */
                            if (this.codon2aa.containsKey(codon)) {
                                String aa = this.codon2aa.get(codon);
                                if (!this.aa2gene2count.containsKey(aa)) {
                                    this.aa2gene2count.put(aa, new HashMap<String, Integer>());
                                }

                                HashMap<String, Integer> gene2count = this.aa2gene2count.get(aa);
                                gene2count.put(gene, gene2count.containsKey(gene) ? gene2count.get(gene) + count : count);
                            } // if corrresponding amino acids exists ... 
                        } // for each count ...
                    } // if counts .length >= 64 
                } // if two parts 
            } // skip the first three lines --
        } // go through the input file line by line --
        br.close();
        return hm;
    } // function .... 

    /**
     * merge two arraylists of double values into a two column RealMatrix
     *
     * @param x
     * @param y
     * @return
     */
    private RealMatrix makeRealMatrix(ArrayList<Double> x, ArrayList<Double> y) {
        double[][] xy = new double[x.size()][2];
        for (int i = 0; i < x.size(); i++) {
            xy[i][0] = x.get(i);
            xy[i][1] = y.get(i);
        }
        return (new Array2DRowRealMatrix(xy));
    }
    
    private ArrayList<String> readFastaForGeneNames(String fastafile) throws FileNotFoundException, IOException {
        ArrayList<String> al = new ArrayList<>();

        // file1 --
        BufferedReader br = new BufferedReader(new FileReader(fastafile));
        String line;
        while ((line = br.readLine()) != null) {
            String[] parts = line.trim().split("\\s+");
            if( parts.length >= 1 && parts[0].startsWith(">") ){
                String genename = parts[0].substring(1);
                al.add(genename);
            }
        }
        br.close();
        return al;
    }

}
