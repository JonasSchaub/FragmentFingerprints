package de.unijena.cheminf.fragment.fingerprint.performanceTest;

import de.unijena.cheminf.fragment.fingerprint.FragmentFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;


public class PerformanceTest {
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    /**
     * Name of file for logging occurred exceptions
     */
    private static final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log.txt";

    /**
     * Name of file for writing results
     */
    private static final String RESULTS_FILE_NAME = "Results";

    /**
     * Name of CSV file for writing time processing time
     */
    private static final String CSV_TIME_FILE_NAME = "CSV_BIT_PROCESS_TIME";

    /**
     * Name of CSV file for the network SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_NETWORK_FILE_NAME = "CSV_COUNT_PROCESS_TIME";

    /**
     * Name of CSV file for the forest SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_FOREST_FILE_NAME = "CSV_Origin_Forest";
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    private ArrayList<ArrayList<String>> deneme = new ArrayList<>();
    private ArrayList<String> fragmentList = new ArrayList<>();
    File dataFile1;
    File dataFile2;
    PrintWriter resultsPrintWriter;
    PrintWriter bitPrintWriter;
    PrintWriter countPrintWriter;

    private ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>();

    public PerformanceTest(String anArgs, String anArgs2) throws IOException {
        /*Set up exception log file*/
        this.workingPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpProcessingTime = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        new File(this.workingPath + "/Results").mkdirs();
        File tmpExceptionsLogFile = new File(this.workingPath + "/Results/"
                + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        PrintWriter tmpExceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        tmpExceptionsPrintWriter.println("#########################################################################");
        tmpExceptionsPrintWriter.println("Processing Time: " + tmpProcessingTime);
        tmpExceptionsPrintWriter.println();
        tmpExceptionsPrintWriter.flush();
        try {
            /*Load SD file*/
             dataFile1 = new File(this.workingPath + anArgs);
             dataFile2 = new File(this.workingPath + anArgs2);
            FileInputStream tmpFile2;
            FileInputStream tmpDBFileInputStream;
            try {
                tmpDBFileInputStream = new FileInputStream(dataFile1);
                tmpFile2 = new FileInputStream(dataFile2);
            } catch (FileNotFoundException | SecurityException anException) {
                throw new IllegalArgumentException("One or more files (name) are invalid: " + anException.getMessage()); // TODO  message
            }
            // results files
            File tmpResultsLogFile = new File(this.workingPath + "/Results/" + PerformanceTest.RESULTS_FILE_NAME + tmpProcessingTime + ".txt");
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            resultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter); // printer
            resultsPrintWriter.println("#########################################################################");
            resultsPrintWriter.println("Processing Time: " + tmpProcessingTime);
            resultsPrintWriter.println("Application initialized. Loading  files named " + anArgs + " and "+ anArgs2+ ".");
           // tmpResultsPrintWriter.flush();
            /*ProcessingTime file*/ // file mit 3 palten
            File tmpCSVTimeFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_TIME_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpCSVTimeFileWriter = new FileWriter(tmpCSVTimeFile, false);
            bitPrintWriter = new PrintWriter(tmpCSVTimeFileWriter);
            //
            File tmpOriginNetworkOriginFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_ORIGIN_NETWORK_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpOriginNetworkOriginFileWriter = new FileWriter(tmpOriginNetworkOriginFile, false);
            countPrintWriter = new PrintWriter(tmpOriginNetworkOriginFileWriter);
           // tmpCSVTimePrintWriter.flush();
            //tmpResultsPrintWriter.println("Hallo");
            // TODO read
            this.generateFingerprints(1,2);
            resultsPrintWriter.flush();
            bitPrintWriter.println();
            bitPrintWriter.flush();
            countPrintWriter.flush();
        } catch (Exception anException) {
            //TODO append log file
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }

    public void read() throws IOException { // TODO exception
        BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader(dataFile1));
        String tmpLine;
        String tmpSeparatorComma = ";";  // TODO separator as argument?
         fragmentList = new ArrayList<>();
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            fragmentList.add(tmpSmilesOfFragments[0]);
        }
        // tmpFragmentSetReader.close();
        fragmentList.remove(0);
        // Read CSV file (itemization tab)
        BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader(dataFile2));
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        tmpMoleculeFragmentsReader.close();
        List<String> tmpSeparateList;
       // ArrayList<ArrayList<String>> deneme = new ArrayList<>();
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
           ArrayList<String> dataForGenerateBitFingerprint = new ArrayList<>();
            for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    dataForGenerateBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                }
            }
            // Illustration the results of the bit arrays for the specified molecules
            // add all molecule maps
            moleculeFragmentList.add(tmpMoleculeFragmentsMap);
            deneme.add(dataForGenerateBitFingerprint);
        }
        /*
        FragmentFingerprinter printer = new FragmentFingerprinter(fragmentList);
        long tmpStartTime = System.currentTimeMillis();
        int i = 0;
        // System.out.println(deneme + "----deneme");
        for (List<String> liste : deneme ) { // deneme sublist für roundNumber
            IBitFingerprint bitFingerprint = printer.getBitFingerprint(liste);
            //  System.out.println(bitFingerprint.cardinality());
            //  i++;
        }
         */
    }
    public void generateFingerprints(long startTime, long endTime) throws IOException {
        this.read();
        FragmentFingerprinter printer = new FragmentFingerprinter(fragmentList);
        resultsPrintWriter.println("Number of molecules: " + deneme.size());
        resultsPrintWriter.println("\n\tGenerate bit fingerprints");
        resultsPrintWriter.println();
        bitPrintWriter.println("Number of fragments: " + fragmentList.size() + " and number of molecules: " + deneme.size());
        bitPrintWriter.println("Number of molecules in process, Bit fingerprint process time in ms");
       // System.out.println(fragmentList + "---fragmentlist");
        for (int i =200000; i <= deneme.size(); i++) { //int i =200; i <= deneme.size(); i+=200
            List<ArrayList<String>> sublist = deneme.subList(0, i);
            startTime = System.currentTimeMillis();
            for (ArrayList<String> a : sublist) {
                IBitFingerprint bit = printer.getBitFingerprint(a);
            }
            endTime = System.currentTimeMillis();
           // System.out.println(endTime-startTime);
            resultsPrintWriter.println("Processing " + sublist.size() + " valid molecules." );
            resultsPrintWriter.println("Bit fingerprint generation took: " + (endTime-startTime) + " ms.");
            bitPrintWriter.println(sublist.size() + "," + (endTime-startTime));
        }
        resultsPrintWriter.println();
        resultsPrintWriter.println("#########################################################################");
        resultsPrintWriter.println("\n\tGenerate count fingerprints");
        resultsPrintWriter.println();
        countPrintWriter.println("Number of molecules in process, Count fingerprint process time in ms.");
        for(int i = 200000; i<= moleculeFragmentList.size(); i++) {  //i+=2000)
            List<HashMap<String,Integer>> sub = moleculeFragmentList.subList(0,i);
            long start = System.currentTimeMillis();
            for(HashMap<String,Integer> map : sub) {
                ICountFingerprint count = printer.getCountFingerprint(map);
            }
            long end = System.currentTimeMillis();
            resultsPrintWriter.println("Processing " + sub.size() + " valid molecules." );
            resultsPrintWriter.println("Count fingerprint generation took: " + (end-start) + " ms.");
            resultsPrintWriter.println();
            countPrintWriter.println(sub.size() + "," + (end-start));
          //  System.out.println(end-start);
        }
       // System.out.println(endTime-startTime + "bit time");



    }


}
