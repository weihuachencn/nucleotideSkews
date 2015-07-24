/*
 * created on April 29, 2014 --
 */
package preferredcodonanalysis;

import com.google.common.base.Joiner;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
import org.apache.commons.cli.ParseException;

/**
 *
 * @author wchen
 * last modified Nov 4, 2014;
 * Weihua Chen ...
 */
public class PreferredCodonAnalysis {

    private final CodonStats cs = new CodonStats(); // codon stats --
    private final ArrayList<String> allCodons = new ArrayList<String>();
    
    public PreferredCodonAnalysis() {
        String codonsString = "TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG";
        for(String codon : codonsString.split("\\s+") ){
            allCodons.add(codon);
        }
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     * @throws org.apache.commons.cli.ParseException
     * @throws java.lang.InterruptedException
     * @throws java.util.concurrent.ExecutionException
     */
    public static void main(String[] args) throws IOException, ParseException, InterruptedException, ExecutionException {
        new PreferredCodonAnalysis().run(args);
    }

    /**
     *
     * @param args
     */
    private void run(String[] args) throws IOException, ParseException, InterruptedException, ExecutionException {
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
        o.addOption("d", "debug", false, "debug mode, only submit 5 jobs; optional");
        o.addOption("n", "threads", true, "number of CPUs to be used; optional; if omitted "
                + "will use available CPUs - 1");
        /**
         * people will be warned anyway ...
         */
        System.err.println("\t!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!\n\tthis program uses a parallel library and it will be VERY memory-consuming\n\tplease try NOT to run it on your desktop / laptop\n\t...so you've been warned\n");
        
        

        /**
         * parse options and check if all required options are there
         */
        CommandLine cmd = parser.parse(o, args);

        if (!cmd.hasOption("c") || !cmd.hasOption("o") ||  !cmd.hasOption("s") || !cmd.hasOption("e")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("this.jar", o);
            System.exit(1);
        }
        
        if( cmd.hasOption("d") ){
            System.err.println("\n\tdebug mode, only the first 20 jobs will be submitted\n\n");
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
         * number of processors to be used; if not set, will use available_CPU -
         * 1 if set but larger than available_CPU, will use all available_CPU if
         * set but available_CPU == 1, will use 1
         */
        int processors = 1;
        int nrOfProcessors = Runtime.getRuntime().availableProcessors();
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

        System.out.println(" ================= your parameters ================= ");
        System.err.println("      use " + processors + " of " + nrOfProcessors + " processor(s)");

        int tasks = 0;
        
        /**
         * paralell stuff
         */
        ExecutorService executor = Executors.newFixedThreadPool(processors);
        ArrayList<Future<String>> futuresList = new ArrayList<Future<String>>();

        for (Map.Entry<String, String> entry : acc2CDSfile.entrySet()) {
            String acc = entry.getKey();
            String fastafile = entry.getValue();
            tasks++;
            PreferredCodonWorkerClass newworker = new PreferredCodonWorkerClass(acc, fastafile, cs, tasks, SeqCountExe, ENCprimeExe, allCodons);
            futuresList.add(executor.submit(newworker));

            if( cmd.hasOption("d") && tasks >= 5 ){
                System.err.println("\n\tstop submitting jobs because you're running in debug mode!!\n\n");
                break; // break the loop !!
            }
        } // for each; to submit jobs --
        executor.shutdown();

        /**
         * prepare output file
         */
        BufferedWriter out;

        out = new BufferedWriter(new FileWriter(outfile));
        String joined = Joiner.on("\t").join( cs.getAAs() );

        out.write(joined);

        out.newLine();

        for (Future<String> result : futuresList) {
            String s = result.get();
            out.write(s); // the results contain new lines already ; Dec 3, 2013 --
        }

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

}
