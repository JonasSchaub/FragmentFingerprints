/*
 * MIT License
 *
 * Copyright (c) 2023 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * A class for testing the performance of the fragment fingerprinter.
 *
 * @author Betuel Sevindik
 */
public class PerformanceTest {
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    /**
     * Name of file for logging occurred exceptions
     */
    private static final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log.txt";
    /**
     * Name of file for writing results
     */
    private static final String RESULTS_FILE_NAME = "Results.txt";
    /**
     * Name of CSV file with the results of the performance test of the bit fingerprints.
     */
    private static final String CSV_BIT_FINGERPRINT_PROCESS_RESULT_FILE_NAME = "CSV_BIT_PROCESS_TIME";
    /**
     * Name of CSV file with the results of the performance test of the count fingerprints.
     */
    private static final String CSV_COUNT_FINGERPRINT_PROCESS_RESULT_FILE_NAME = "CSV_COUNT_PROCESS_TIME";
    /**
     * Name of the CSV file with the results of the generated bit arrays.
     */
    private static final String BIT_ARRAY_RESULT_FILE_NAME = "BIT_ARRAY_TIME";
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Initial capacity of the lists in which the data for generating the fingerprints are stored.
     */
    private final int INITIAL_CAPACITY_VALUE = 200000;
    /**
     * Separator for separating the lines in the input files
     */
    private final String LINE_SEPARATOR_COMMA = ",";
    /**
     * Separator for separating the lines in the input files.
     */
    private final String LINE_SEPARATOR_SEMICOLON = ";";
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private  class variables">
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    /**
     * List which is contains all the data of the molecules needed to generate the bit fingerprint.
     */
    private ArrayList<ArrayList<String>> moleculeListForBitFingerprint = new ArrayList<>(this.INITIAL_CAPACITY_VALUE);
    /**
     * List that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private ArrayList<String> fragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE);
    /**
     * Stores the name/ID of the molecules from moleculeFile
     */
    private List<String> listOfMoleculeNames;
    /**
     * CSV file containing the resulting fragments and their frequencies for each given SMILES (molecule).
     */
    private File moleculeFile;
    /**
     * CSV file containing the resulting fragments.
     */
    private File fragmentFile;
    /**
     * PrintWriter for logging the results of fingerprint generation.
     */
    private PrintWriter resultsPrintWriter;
    /**
     * PrintWriter for logging the bit fingerprint generation.
     */
    private PrintWriter bitPrintWriter;
    /**
     * PrintWriter for logging the count fingerprint generation.
     */
    private PrintWriter countPrintWriter;
    /**
     * PrintWriter for logging exceptions.
     */
    private PrintWriter exceptionsPrintWriter;
    /**
     * PrintWriter to log the generated bit arrays.
     */
    private PrintWriter bitArrayPrintWriter;
    /**
     * list stores each molecule as HashMap. The HashMap assigns
     * the unique SMILES (fragments) of the molecule to the frequency of the molecule,
     * which is read out from the molecule file.
     */
    private ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE);
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor.
     * All files listed below as parameters must be located in the same directory as the application JAR file.
     * In the constructor all necessary files for storing the performance test are created.
     * There are 3 different files to document the result of the performance snapshot.
     * The results during the generation of the bit and count fingerprint are stored once
     * individually and once together.
     * Subsequently, the fingerprints are created, and the results are documented in the corresponding files.
     *
     * @param anArgs the command line arguments, anArgs[0] and anArgs[1] must be the names of the text files to load.
     *               The text files must be located in the same directory as the application's JAR file. The command
     *               line anArgs[3] must be the name of the path, where the performance test results are
     *               stored. The command line anArgs[2] must be the number of the molecules in process.
     * @throws IOException is thrown if the constructor is unable to open a text file for logging occurred exceptions.
     */
    public PerformanceTest(String[] anArgs) throws IOException {
        if(anArgs.length !=4) {
            throw new IllegalArgumentException("Four arguments (three file names and the number of the molecules in process) are required.");
        }
        String tmpFilePath =  (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        /*Set up exception log file*/
        this.workingPath = (new File(anArgs[2]).getAbsoluteFile().getAbsolutePath()) + File.separator;
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
            this.moleculeFile = new File(tmpFilePath + anArgs[0]);
            this.fragmentFile = new File(tmpFilePath + anArgs[1]);
            FileInputStream tmpFragmentFileInputStream;
            FileInputStream tmpMoleculeFileInputStream;
            try {
                tmpMoleculeFileInputStream = new FileInputStream(this.moleculeFile);
                tmpFragmentFileInputStream = new FileInputStream(this.fragmentFile);
            } catch (FileNotFoundException | SecurityException anException) {
                this.appendToLogfile(anException);
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
            this.resultsPrintWriter.println("Application initialized. Loading  files named " + anArgs[0] + " and "+ anArgs[1] + ".");
            this.resultsPrintWriter.println();
            // bit fingerprints process file
            File tmpBitFingerprintsResultFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_BIT_FINGERPRINT_PROCESS_RESULT_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpBitFingerprintsResultWriter = new FileWriter(tmpBitFingerprintsResultFile, false);
            this.bitPrintWriter = new PrintWriter(tmpBitFingerprintsResultWriter);
            // count fingerprints process file
            File tmpCountFingerprintsResultFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_COUNT_FINGERPRINT_PROCESS_RESULT_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpCountFingerprintResultWriter = new FileWriter(tmpCountFingerprintsResultFile, false);
            this.countPrintWriter = new PrintWriter(tmpCountFingerprintResultWriter);
            // bit array file
            File tmpBitArrayFingerprintResultFile = new File(this.workingPath + "/Results/" + PerformanceTest.BIT_ARRAY_RESULT_FILE_NAME + tmpProcessingTime+ ".txt");
            FileWriter tmpBitArrayFingerprintResultWriter = new FileWriter(tmpBitArrayFingerprintResultFile, false);
            this.bitArrayPrintWriter = new PrintWriter(tmpBitArrayFingerprintResultWriter);
            // read in CSV files that contains fragments
            try {
                this.importDataFromTextFile();
            } catch (IOException anException) {
                this.exceptionsPrintWriter.println("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in.");
                this.exceptionsPrintWriter.flush();
                this.appendToLogfile(anException);
                throw new Exception("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in. ");
            }
            // generate bit and count fingerprints
            this.generateFingerprints(Integer.parseInt(anArgs[3]),1,2); // TODO delete arguments
            this.resultsPrintWriter.flush();
            this.bitPrintWriter.println();
            this.bitPrintWriter.flush();
            this.countPrintWriter.flush();
            this.bitArrayPrintWriter.flush();
            this.exceptionsPrintWriter.flush();
        } catch (Exception anException) {
            this.appendToLogfile(anException);
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * The two text files (fragment and molecule file)
     * are read in and the data provided so that in the next step the
     * fingerprints can be calculated based on the information in these text files.
     * For the method to work properly, it is important that the two files are in a specific form predicted
     * for the method. Thus, each line in the fragment file corresponds to a fragment. Components
     * within a line must be separated by
     * a semicolon. The fragment represented by a unique SMILES must be at the very front (first position) of the line.
     * All other information within the line is ignored.
     *
     * The same applies to the molecule file. Here, too, a molecule is defined in each line. However, each line
     * starts with the name/ID of the molecule and the SMILES encoding corresponding to the molecule structure.
     * Then, the fragments of the molecule and the corresponding frequency of the fragment in the molecule follow
     * alternately in the line. Again, all components are separated by a semicolon.
     *
     * @throws IOException is thrown if an error occurs when reading in the two text files.
     */
    private void importDataFromTextFile() throws IOException {
        BufferedReader tmpFragmentSetReader;
        BufferedReader tmpMoleculeFragmentsReader;
        try {
            tmpFragmentSetReader =  new BufferedReader(new FileReader(this.moleculeFile));
            tmpMoleculeFragmentsReader = new BufferedReader(new FileReader(this.fragmentFile));
        } catch(IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("File is not readable");
        }
        //Read CSV file (fragments file) to obtain fragments used to create fingerprints
        String tmpLine;
        this.fragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE);
        try {
            while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                String[] tmpSmilesOfFragments = tmpLine.split(this.LINE_SEPARATOR_SEMICOLON);
                this.fragmentList.add(tmpSmilesOfFragments[0]);
            }
            this.fragmentList.remove(0);
        } catch(IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid fragment file. At least one line is not readable.");
        }
        // Read CSV file (molecules file)
        List<List<String>> tmpList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE);
        String tmpMoleculeLine;
        try {
            while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
                String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(this.LINE_SEPARATOR_SEMICOLON);
                tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
            }
        } catch (IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid molecule file. At least one line is not readable");
        }
        tmpMoleculeFragmentsReader.close();
        List<String> tmpSeparateList;
        this.listOfMoleculeNames = new ArrayList<>();
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            this.listOfMoleculeNames.add(tmpSeparateList.get(0));
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            ArrayList<String> dataForGenerateBitFingerprint = new ArrayList<>();
            try {
                for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                    if (i % 2 == 0) {
                        tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                        dataForGenerateBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                    }
                }
            } catch (IndexOutOfBoundsException anException) {
                this.appendToLogfile(anException);
                throw new IndexOutOfBoundsException("The line has not the right length");
            }
            this.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
            this.moleculeListForBitFingerprint.add(dataForGenerateBitFingerprint);
        }
    }
    //
    /**
     * Starts the generation of the fingerprints. And writes the results into the corresponding text file.
     *
     * @param aNumberOfMoleculesInProcess specifies how many fingerprints are to be generated.
     * @param aEndTime System time in ms.
     * @param aStartTime System time in ms.
     * @throws Exception is thrown if an error occurs when generating the fingerprints.
     */
    private void generateFingerprints(int aNumberOfMoleculesInProcess, long aEndTime, long aStartTime) throws Exception {
        FragmentFingerprinter tmpFragmentFingerprinter = new FragmentFingerprinter(this.fragmentList);
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("Number of molecules: " + this.moleculeListForBitFingerprint.size());
        this.resultsPrintWriter.println("Number of fragment: " + this.moleculeListForBitFingerprint.size());
        this.resultsPrintWriter.println("Succeeded loading of molecules and fragments");
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("\n\tGenerate bit fingerprints");
        this.resultsPrintWriter.println();
        this.bitPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListForBitFingerprint.size());
        this.bitPrintWriter.println("Number of molecules in process, Bit fingerprint process time in ms");
        this.bitArrayPrintWriter.println("Molecule name/ID, bit array");
        for (int i = aNumberOfMoleculesInProcess; i <= this.moleculeListForBitFingerprint.size(); i+= aNumberOfMoleculesInProcess) {
            List<ArrayList<String>> tmpNumberOfMoleculesInProcess = this.moleculeListForBitFingerprint.subList(0, i);
            aStartTime = System.currentTimeMillis();
            for (ArrayList<String> tmpMolecule : tmpNumberOfMoleculesInProcess) {
                try {
                    IBitFingerprint bit = tmpFragmentFingerprinter.getBitFingerprint(tmpMolecule);
                } catch (Exception anException) {
                    this.exceptionsPrintWriter.println("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                    this.appendToLogfile(anException);
                    throw new Exception("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                }
            }
            aEndTime = System.currentTimeMillis();
            this.resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcess.size() + " valid molecules.");
            this.resultsPrintWriter.println("Bit fingerprint generation took: " + (aEndTime - aStartTime) + " ms.");
            this.bitPrintWriter.println(tmpNumberOfMoleculesInProcess.size() + "," + (aEndTime - aStartTime));
        }
        int[] tmpBitArray = new int[this.fragmentList.size()];
        for(Iterator tmpMoleculeNameIterator = this.listOfMoleculeNames.iterator(), tmpMolecule = this.moleculeListForBitFingerprint.iterator(); tmpMoleculeNameIterator.hasNext() && tmpMolecule.hasNext();) { //ArrayList<String> tmpMolecule : this.moleculeListForBitFingerprint
            List<String> tmpListOfMolecule = (List<String>) tmpMolecule.next();
            String tmpMoleculeNameOrID = (String) tmpMoleculeNameIterator.next();
            tmpBitArray = tmpFragmentFingerprinter.getBitArray(tmpListOfMolecule);
            this.bitArrayPrintWriter.println(tmpMoleculeNameOrID + ":" + java.util.Arrays.toString(tmpBitArray));
        }
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("#########################################################################");
        this.resultsPrintWriter.println("\n\tGenerate count fingerprints");
        this.resultsPrintWriter.println();
        this.countPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListForBitFingerprint.size());
        this.countPrintWriter.println("Number of molecules in process, Count fingerprint process time in ms.");
        for(int i = aNumberOfMoleculesInProcess; i<= this.moleculeFragmentList.size(); i+= aNumberOfMoleculesInProcess) {  //i+=2000)
            List<HashMap<String, Integer>> tmpNumberOfMoleculesInProcessList = this.moleculeFragmentList.subList(0, i);
            long start = System.currentTimeMillis();
            for (HashMap<String, Integer> tmpUniqueSmilesToFrequencyMap : tmpNumberOfMoleculesInProcessList) {
                try {
                    ICountFingerprint count = tmpFragmentFingerprinter.getCountFingerprint(tmpUniqueSmilesToFrequencyMap);
                } catch (Exception anException) {
                    this.exceptionsPrintWriter.println("Count fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                    this.appendToLogfile(anException);
                    throw new Exception("Count fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                }
            }
            long end = System.currentTimeMillis();
            this.resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcessList.size() + " valid molecules.");
            this.resultsPrintWriter.println("Count fingerprint generation took: " + (end - start) + " ms.");
            this.resultsPrintWriter.println();
            this.countPrintWriter.println(tmpNumberOfMoleculesInProcessList.size() + "," + (end - start));
        }
    }
    //
    /**
     * Appends the given exception's stack trace to a log file.
     *
     * @param anException the exception to log
     */
    private void appendToLogfile(Exception anException) {
        if (anException == null) {
            return;
        }
        this.exceptionsPrintWriter = null;
        PrintWriter m = null;
        try {
            FileWriter tmpFileWriter = new FileWriter(this.workingPath
                    + "/Results/" + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME, true);
            this.exceptionsPrintWriter = new PrintWriter(tmpFileWriter);
            StringWriter tmpStringWriter = new StringWriter();
            anException.printStackTrace(new PrintWriter(tmpStringWriter));
            String tmpStackTrace = tmpStringWriter.toString();
            this.exceptionsPrintWriter.println(tmpStackTrace);
        } catch (IOException anIOException) {
            anIOException.printStackTrace(System.err);
        } finally {
            if (this.exceptionsPrintWriter != null) {
                this.exceptionsPrintWriter.println();
                this.exceptionsPrintWriter.flush();
                this.exceptionsPrintWriter.close();
            }
        }
    }
    //</editor-fold>
}
