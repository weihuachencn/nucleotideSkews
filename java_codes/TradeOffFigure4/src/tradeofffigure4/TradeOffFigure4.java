/*
 */
package tradeofffigure4;

import com.google.common.base.Joiner;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

/**
 * created : Oct 17, 2013 last modified : Apr 5, 2014 last modified : May 22,
 * 2014; added AA costs under other (unusual) conditions --
 *
 * June 15, 2014; added NSP -- June 16, 2014: added options for gcstart, gcend
 * 
 * Sep 21, 2014; added AT skews only or GC skews only 
 *
 * @author wchen
 */
public class TradeOffFigure4 {

    private int tasks = 0, executedTasks = 0;
    private BufferedWriter out;
    private long begTest, endTest;
    // codon stats table --
    private HashMap<String, HashMap<Integer, Integer>> hmCodonStates;
    private ArrayList<HashMap<String, Float>> arAA2Cost;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, ParseException, InterruptedException, ExecutionException, org.apache.commons.cli.ParseException {
        new TradeOffFigure4().run(args);
    }

    public TradeOffFigure4() {
        begTest = new java.util.Date().getTime();
        hmCodonStates = this.codonStats();
        arAA2Cost = this.getAA2Cost();
    }

    /**
     * run
     *
     * @param args
     */
    private void run(String[] args) throws IOException, ParseException, InterruptedException, ExecutionException, org.apache.commons.cli.ParseException {

        int nrOfProcessors = Runtime.getRuntime().availableProcessors();

        System.out.println("--> in total " + nrOfProcessors + " processers available on your machine");
        System.out.println();

        /**
         * command line parser using CLI.commons
         */
        // create the command line parser
        CommandLineParser parser = new BasicParser();

        // create the Options
        Options o = new Options();

        o.addOption("o", "output-file", true, "output file; mandatory");
        o.addOption("n", "threads", true, "number of CPUs to be used; if omitted "
                + "will use available CPUs - 1");
        o.addOption("d", "debug", false, "debug mode, only submit 5 jobs");

        o.addOption("a", "gcstart", true, "GC to start with; integer; optional, default = 10");
        o.addOption("b", "gcend", true, "GC to end with; integer; optional, default = 90");

        o.addOption("x", "biasstart", true, "bias to start with; integer; optional, default = 10");
        o.addOption("y", "biasend", true, "bias to end with; integer; optional, default = 90");
        
        // at or gc skews only 
        o.addOption("s", "skewsonly", true, "at or gc skews only; string; optional, default = all; acceptable options are : all|atonly|gconly");

        System.out.println("--> there will be about 13,600 tasks to be done; \n"
                + "--> it will take about a few hours using 10-CUP cores; \n"
                + "--> and it takes a huge amount of memory (>=20G); \n"
                + "--> have fun!! \n\n");

        /**
         * parse options and check if all required options are there
         */
        CommandLine cmd = parser.parse(o, args);

        if (!cmd.hasOption("o")) {

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("java -Xmx6G -jar \"/directory/to/this.jar\"", o);
            System.exit(1);
        }
        
        String skews = "all";
        if( cmd.hasOption("s") ){
            skews = cmd.getOptionValue("s");
            if( !skews.equalsIgnoreCase("all") &&  !skews.equalsIgnoreCase("atonly") && !skews.equalsIgnoreCase("gconly")){
                System.err.println("invalid option value for 's':" + skews);
                System.err.println("valid option value should be one of the following: " );
                System.err.println("  all, atonly, gconly" );
                System.err.println("now quiting ..." );
                System.exit(1);
            }
        }
        
        String outfile = cmd.getOptionValue("o");

        /**
         * number of processors to be used; if not set, will use available_CPU -
         * 1 if set but larger than available_CPU, will use all available_CPU if
         * set but available_CPU == 1, will use 1
         */
        int processors = 1;
        if (cmd.hasOption("n")) {
            try {
                processors = Integer.parseInt(cmd.getOptionValue("n"));

                if (processors > nrOfProcessors) {
                    processors = nrOfProcessors;
                }
            } catch (NumberFormatException e) {
                System.err.println("illegal option -n");
            }
        } else {
            processors = nrOfProcessors > 1 ? nrOfProcessors - 1 : 1;
        }

        /**
         * start and end GCs --
         */
        int startGC = 10, endGC = 90;
        if (cmd.hasOption("a")) {
            try {
                startGC = Integer.parseInt(cmd.getOptionValue("a"));

                if (startGC < 10) {
                    startGC = 10;
                }
            } catch (NumberFormatException e) {
                System.err.println("illegal option -a");
            }
        }

        if (cmd.hasOption("b")) {
            try {
                endGC = Integer.parseInt(cmd.getOptionValue("b"));

                if (endGC > 90) {
                    endGC = 90;
                }
            } catch (NumberFormatException e) {
                System.err.println("illegal option -b");
            }
        }

        if (startGC >= endGC) {
            startGC = 10;
            endGC = 90;
        }

        /**
         * start and end bias --
         */
        int startbias = 10, endbias = 90;
        if (cmd.hasOption("x")) {
            try {
                startbias = Integer.parseInt(cmd.getOptionValue("x"));

                if (startbias < 10) {
                    startbias = 10;
                }
            } catch (NumberFormatException e) {
                System.err.println("illegal option -x");
            }
        }

        if (cmd.hasOption("y")) {
            try {
                endbias = Integer.parseInt(cmd.getOptionValue("y"));

                if (endbias > 90) {
                    endbias = 90;
                }
            } catch (NumberFormatException e) {
                System.err.println("illegal option -y");
            }
        }

        if (startbias >= endbias) {
            startbias = 10;
            endbias = 90;
        }

        System.out.println(" ================= your parameters ================= ");
        System.err.println("      use " + processors + " of " + nrOfProcessors + " processor(s)");
        System.err.println("      GC range: " + startGC + " ~ " + endGC);
        System.err.println("      bias range: " + startbias + " ~ " + endbias);

        /**
         * paralelle stuff
         */
        ExecutorService executor = Executors.newFixedThreadPool(processors);
        ArrayList<Future<String>> futuresList = new ArrayList<Future<String>>();

        for (float gc = startGC; gc <= endGC; gc += 0.1) {
            for (int bias = startbias; bias <= endbias; bias += 5) {
                tasks++;
                TradeoffWorkerClass newSimTask = new TradeoffWorkerClass(gc, bias, hmCodonStates, arAA2Cost, tasks, skews);
                futuresList.add(executor.submit(newSimTask));

                // debug mode
                if (cmd.hasOption("d") && tasks >= 5) {
                    System.err.println("\n\tstop submitting jobs because you're running in debug mode!!\n\n");
                    break; // break the loop !!
                }
            }

            // if debug mode !!
            if (cmd.hasOption("d") && tasks >= 5) {
                System.err.println("\n\tstop submitting jobs because you're running in debug mode!!\n\n");
                break; // break the loop !!
            }
        }
        executor.shutdown();

        /**
         * prepare output file
         */
        out = new BufferedWriter(new FileWriter(outfile));
        String joined = Joiner.on("\t").join(Arrays.asList("GC", "overallBiases", "overallSkews", "AT1", "GC1", "AT2", "GC2", "AT4", "GC4", "aaAvgCost",
                "aaAvgCost_Aer_Het", "aaAvgCost_An_Het", "aaAvgCost_Aer_Phot", "aaAvgCost_An_Phot", "aaAVGNSP", "Carbon"));
        out.write(joined);
        out.newLine();

        for (Future<String> result : futuresList) {
            String s = result.get();
            out.write(s); // the results contain new lines already ; Dec 3, 2013 --
        }
        out.close();

        // all jobs done!! 
    }// run

    /**
     * @return $hash{ $codon }{ 1,2,3 } = 1,2,3,4
     */
    private HashMap<String, HashMap<Integer, Integer>> codonStats() {

        HashMap<String, HashMap<Integer, Integer>> hmResult = new HashMap<>();
        HashMap<String, String> hmCodon2AA = new HashMap<>();

        // codons
        char[] aB1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".toCharArray();
        char[] aB2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".toCharArray();
        char[] aB3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".toCharArray();
        char[] aA = "FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".toCharArray();

        // get codon to aa hashmap; stop codons are skipped --
        for (int i = 0; i < aB1.length; i++) {
            String aa = aA[i] + "";
            StringBuilder codon = new StringBuilder();
            codon.append(aB1[i]).append(aB2[i]).append(aB3[i]);
            if (!aa.equalsIgnoreCase("X") && codon.length() == 3) {
                String codonStr = codon.toString();
                hmCodon2AA.put(codonStr, aa);
            }
        }// 

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
                    String newaa = hmCodon2AA.containsKey(ncodon) ? hmCodon2AA.get(ncodon) : "stop";
                    if (newaa.equalsIgnoreCase(aa)) {
                        dindex++;
                    }
                } // end of for each nucleotide

                // -- --
                if (!hmResult.containsKey(codon)) {
                    hmResult.put(codon, new HashMap<Integer, Integer>());
                }
                hmResult.get(codon).put(p, dindex);
            }// iterate each codon position
        }// iterate each codon
        return hmResult;
    } // end of codonStats

    private ArrayList<HashMap<String, Float>> getAA2Cost() {
        ArrayList<HashMap<String, Float>> arhash = new ArrayList<>();

        /**
         * will calcualte all at the same time --
         */
        // 0, the pnas paper --
        HashMap<String, Float> zero = new HashMap<>();
        zero.put("A", 11.7f);
        zero.put("G", 11.7f);
        zero.put("S", 11.7f);
        zero.put("D", 12.7f);
        zero.put("N", 14.7f);
        zero.put("E", 15.3f);
        zero.put("Q", 16.3f);
        zero.put("T", 18.7f);
        zero.put("P", 20.3f);
        zero.put("V", 23.3f);
        zero.put("C", 24.7f);
        zero.put("L", 27.3f);
        zero.put("R", 27.3f);
        zero.put("K", 30.3f);
        zero.put("I", 32.3f);
        zero.put("M", 34.3f);
        zero.put("H", 38.3f);
        zero.put("Y", 50f);
        zero.put("F", 52f);
        zero.put("W", 74.3f);

        arhash.add(zero);

        // 1, the jbe paper aerobic heterotroph (1) (PMID: 22538926)
        HashMap<String, Float> one = new HashMap<>();
        one.put("A", 14.5f);
        one.put("R", 20.5f);
        one.put("N", 18.5f);
        one.put("D", 15.5f);
        one.put("C", 26.5f);
        one.put("Q", 10.5f);
        one.put("E", 9.5f);
        one.put("G", 14.5f);
        one.put("H", 29f);
        one.put("I", 38f);
        one.put("L", 37f);
        one.put("K", 36f);
        one.put("M", 36.5f);
        one.put("F", 61f);
        one.put("P", 14.5f);
        one.put("S", 14.5f);
        one.put("T", 21.5f);
        one.put("W", 75.5f);
        one.put("Y", 59f);
        one.put("V", 29f);

        arhash.add(one);

        // 2,  anaerobic heterotroph (2) (PMID: 22538926)
        HashMap<String, Float> two = new HashMap<>();
        two.put("A", 2f);
        two.put("R", 13f);
        two.put("N", 6f);
        two.put("D", 3f);
        two.put("C", 13f);
        two.put("Q", 3f);
        two.put("E", 2f);
        two.put("G", 1f);
        two.put("H", 5f);
        two.put("I", 14f);
        two.put("L", 4f);
        two.put("K", 12f);
        two.put("M", 24f);
        two.put("F", 10f);
        two.put("P", 7f);
        two.put("S", 1f);
        two.put("T", 9f);
        two.put("W", 14f);
        two.put("Y", 8f);
        two.put("V", 4f);
        arhash.add(two);

        // Photaerobic phototroph (3) (PMID: 22538926)");
        HashMap<String, Float> three = new HashMap<>();
        three.put("A", 14.5f);
        three.put("R", 20.5f);
        three.put("N", 18.5f);
        three.put("D", 15.5f);
        three.put("C", 26.5f);
        three.put("Q", 10.5f);
        three.put("E", 9.5f);
        three.put("G", 14.5f);
        three.put("H", 31f);
        three.put("I", 38f);
        three.put("L", 37f);
        three.put("K", 36f);
        three.put("M", 36.5f);
        three.put("F", 63f);
        three.put("P", 14.5f);
        three.put("S", 14.5f);
        three.put("T", 21.5f);
        three.put("W", 77.5f);
        three.put("Y", 61f);
        three.put("V", 29f);

        arhash.add(three);

        // anaerobic phototroph (4) (PMID: 22538926)");
        HashMap<String, Float> four = new HashMap<>();
        four.put("A", 2f);
        four.put("R", 13f);
        four.put("N", 6f);
        four.put("D", 3f);
        four.put("C", 13f);
        four.put("Q", 3f);
        four.put("E", 2f);
        four.put("G", 1f);
        four.put("H", 7f);
        four.put("I", 14f);
        four.put("L", 4f);
        four.put("K", 12f);
        four.put("M", 24f);
        four.put("F", 14f);
        four.put("P", 7f);
        four.put("S", 1f);
        four.put("T", 9f);
        four.put("W", 16f);
        four.put("Y", 10f);
        four.put("V", 4f);

        arhash.add(four);

        // NSP; June 15, 2014 --
        HashMap<String, Float> five = new HashMap<>();
        five.put("A", 1f);
        five.put("R", 4f);
        five.put("N", 2f);
        five.put("D", 1f);
        five.put("C", 2f);
        five.put("E", 2f);
        five.put("Q", 1f);
        five.put("G", 1f);
        five.put("H", 3f);
        five.put("I", 1f);
        five.put("L", 1f);
        five.put("K", 2f);
        five.put("M", 2f);
        five.put("F", 1f);
        five.put("P", 1f);
        five.put("S", 1f);
        five.put("T", 1f);
        five.put("W", 2f);
        five.put("Y", 1f);
        five.put("V", 1f);

        arhash.add(five);

        // Number of Carbons per amino acid 
        HashMap<String, Float> six = new HashMap<>();
        six.put("A", 3f);
        six.put("R", 6f);
        six.put("N", 4f);
        six.put("D", 4f);
        six.put("C", 3f);
        six.put("E", 5f);
        six.put("Q", 5f);
        six.put("G", 2f);
        six.put("H", 6f);
        six.put("I", 6f);
        six.put("L", 6f);
        six.put("K", 6f);
        six.put("M", 5f);
        six.put("F", 8f);
        six.put("P", 5f);
        six.put("S", 3f);
        six.put("T", 4f);
        six.put("W", 11f);
        six.put("Y", 9f);
        six.put("V", 5f);

        arhash.add(six);

        // return 
        return arhash;
    }
}
