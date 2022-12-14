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
import java.io.StringWriter;
import java.text.NumberFormat;
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
    private ArrayList<ArrayList<String>> moleculeListforBitFingerprint = new ArrayList<>();
    private ArrayList<String> fragmentList = new ArrayList<>();
    File moleculeFile;
    File fragmentFile;
    PrintWriter resultsPrintWriter;
    PrintWriter bitPrintWriter;
    PrintWriter countPrintWriter;
    PrintWriter exceptionsPrintWriter;

    private ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>();

    public PerformanceTest(String anArgs, String anArgs2, String anArgs3) throws IOException {
        /*Set up exception log file*/
        this.workingPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpProcessingTime = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        new File(this.workingPath + "/Results").mkdirs();
        File tmpExceptionsLogFile = new File(this.workingPath + "/Results/"
                + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        this.exceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        this.exceptionsPrintWriter.println("#########################################################################");
        this.exceptionsPrintWriter.println("Processing Time: " + tmpProcessingTime);
        this.exceptionsPrintWriter.println();
        this.exceptionsPrintWriter.flush();
        try {
            // Load fragment files
            this.moleculeFile = new File(this.workingPath + anArgs);
            this.fragmentFile = new File(this.workingPath + anArgs2);
            FileInputStream tmpFragmentFileInputStream;
            FileInputStream tmpMoleculeFileInputStream;
            try {
                tmpMoleculeFileInputStream = new FileInputStream(moleculeFile);
                tmpFragmentFileInputStream = new FileInputStream(fragmentFile);
            } catch (FileNotFoundException | SecurityException anException) {
                throw new IllegalArgumentException("One or more files (name) are invalid: " + anException.getMessage());
            }
            // results files
            File tmpResultsLogFile = new File(this.workingPath + "/Results/" + PerformanceTest.RESULTS_FILE_NAME + tmpProcessingTime + ".txt");
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            this.resultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            this.resultsPrintWriter.println("#########################################################################");
            this.resultsPrintWriter.println();
            this.resultsPrintWriter.println("Processing Time: " + tmpProcessingTime);
            this.resultsPrintWriter.println();
            this.resultsPrintWriter.println("Application initialized. Loading  files named " + anArgs + " and "+ anArgs2+ ".");
            this.resultsPrintWriter.println();
            // memory usages
            this.resultsPrintWriter.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " / " + Runtime.getRuntime().totalMemory()/(1024*1024) + " Mb");
            // bit fingerprints process file
            File tmpCSVTimeFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_TIME_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpCSVTimeFileWriter = new FileWriter(tmpCSVTimeFile, false);
            this.bitPrintWriter = new PrintWriter(tmpCSVTimeFileWriter);
            // count fingerprints process file
            File tmpOriginNetworkOriginFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_ORIGIN_NETWORK_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpOriginNetworkOriginFileWriter = new FileWriter(tmpOriginNetworkOriginFile, false);
            this.countPrintWriter = new PrintWriter(tmpOriginNetworkOriginFileWriter);
            // read in CSV files that contains fragments
            this.generateFingerprints(Integer.parseInt(anArgs3),1,2);
            this.resultsPrintWriter.flush();
            this.bitPrintWriter.println();
            this.bitPrintWriter.flush();
            this.countPrintWriter.flush();
            this.exceptionsPrintWriter.flush();
        } catch (Exception anException) {
            this.appendToLogfile(anException);
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }

    public void LoadData() throws IOException {
        BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader(this.moleculeFile));
        String tmpLine;
        String tmpSeparatorComma = ";";  // TODO separator as argument?
        this.fragmentList = new ArrayList<>();
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            this.fragmentList.add(tmpSmilesOfFragments[0]);
        }
        this.fragmentList.remove(0);
        // Read CSV file (itemization tab)
        BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader(this.fragmentFile));
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        tmpMoleculeFragmentsReader.close();
        List<String> tmpSeparateList;
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
            this.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
            this.moleculeListforBitFingerprint.add(dataForGenerateBitFingerprint);
        }
    }
    public void generateFingerprints(int numberOfMoleculesInProcess, long endTime, long startTime) throws Exception {
        try {
            this.LoadData();
        } catch (Exception anException) {
            this.exceptionsPrintWriter.println("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in.");
            this.exceptionsPrintWriter.flush();
            this.appendToLogfile(anException);
            throw new Exception("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in. ");
        }
        FragmentFingerprinter printer = new FragmentFingerprinter(this.fragmentList);
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("Number of molecules: " + this.moleculeListforBitFingerprint.size());
        this.resultsPrintWriter.println("Number of fragment: " + this.moleculeListforBitFingerprint.size());
        this.resultsPrintWriter.println("Succeeded loading of molecules and fragments");
        this.resultsPrintWriter.println();
        NumberFormat tmpNumberFormat = NumberFormat.getNumberInstance();
        this.resultsPrintWriter.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " / " + Runtime.getRuntime().totalMemory()/(1024*1024) + " Mb");
       // System.out.println("Memory usage: " + tmpNumberFormat.format((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " MB" + " / " + tmpNumberFormat.format(Runtime.getRuntime().totalMemory()/(1024* 1024) + " MB")));
        this.resultsPrintWriter.println("\n\tGenerate bit fingerprints");
        this.resultsPrintWriter.println();
        this.bitPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListforBitFingerprint.size());
        this.bitPrintWriter.println("Number of molecules in process, Bit fingerprint process time in ms");
        for (int i = numberOfMoleculesInProcess; i <= this.moleculeListforBitFingerprint.size(); i+=numberOfMoleculesInProcess) { // +=numberOfMoleculesInProcess
                List<ArrayList<String>> tmpNumberOfMoleculesInProcess = this.moleculeListforBitFingerprint.subList(0, i);
                try {
                 startTime = System.currentTimeMillis();
                    for (ArrayList<String> tmpMolecule : tmpNumberOfMoleculesInProcess) {
                        IBitFingerprint bit = printer.getBitFingerprint(tmpMolecule);
                    }
                     endTime = System.currentTimeMillis();
                    /*
                    this.resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcess.size() + " valid molecules.");
                    this.resultsPrintWriter.println("Bit fingerprint generation took: " + (endTime - startTime) + " ms.");
                    this.bitPrintWriter.println(tmpNumberOfMoleculesInProcess.size() + "," + (endTime - startTime));

                     */
                } catch (Exception anException) {
                    this.exceptionsPrintWriter.println("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list or molecule list.");
                    this.appendToLogfile(anException);
                    throw new Exception("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list or molecule list.");
                }
                resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcess.size() + " valid molecules.");
                resultsPrintWriter.println("Bit fingerprint generation took: " + (endTime - startTime) + " ms.");
                bitPrintWriter.println(tmpNumberOfMoleculesInProcess.size() + "," + (endTime - startTime));
        }
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " / " + Runtime.getRuntime().totalMemory()/(1024*1024) + " Mb");
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("#########################################################################");
        this.resultsPrintWriter.println("\n\tGenerate count fingerprints");
        this.resultsPrintWriter.println();
        this.countPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListforBitFingerprint.size());
        this.countPrintWriter.println("Number of molecules in process, Count fingerprint process time in ms.");
        for(int i = numberOfMoleculesInProcess; i<= this.moleculeFragmentList.size(); i+=numberOfMoleculesInProcess) {  //i+=2000)
                List<HashMap<String, Integer>> sub = this.moleculeFragmentList.subList(0, i);
                try {
                    long start = System.currentTimeMillis();
                    for (HashMap<String, Integer> map : sub) {
                        ICountFingerprint count = printer.getCountFingerprint(map);
                    }
                    long end = System.currentTimeMillis();
                    this.resultsPrintWriter.println("Processing " + sub.size() + " valid molecules.");
                    this.resultsPrintWriter.println("Count fingerprint generation took: " + (end - start) + " ms.");
                    this.resultsPrintWriter.println();
                    this.countPrintWriter.println(sub.size() + "," + (end - start));
                } catch (Exception anException) {
                    this.exceptionsPrintWriter.println("Count fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                    this.appendToLogfile(anException);
                    throw new Exception("Count fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                }
            }
        this.resultsPrintWriter.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " / " + Runtime.getRuntime().totalMemory()/(1024*1024) + " Mb");
        }
    private void appendToLogfile(Exception anException) {
        if (anException == null) {
            return;
        }
        PrintWriter tmpPrintWriter = null;
        try {
            FileWriter tmpFileWriter = new FileWriter(this.workingPath
                    + "/Results/" + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME, true);
            tmpPrintWriter = new PrintWriter(tmpFileWriter);
            StringWriter tmpStringWriter = new StringWriter();
            anException.printStackTrace(new PrintWriter(tmpStringWriter));
            String tmpStackTrace = tmpStringWriter.toString();
            tmpPrintWriter.println(tmpStackTrace);
        } catch (IOException anIOException) {
            anIOException.printStackTrace(System.err);
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.println();
                tmpPrintWriter.flush();
                tmpPrintWriter.close();
            }
        }
    }
}
