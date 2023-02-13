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
import java.util.List;

/**
 * A test class for testing the performance of the FragmentFingerprinter methods.
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
    private static final String RESULTS_FILE_NAME = "Results";
    /**
     * Name of CSV file with the results of the performance test of the bit fingerprints.
     */
    private static final String CSV_BIT_FINGERPRINT_PROCESS_RESULT = "CSV_BIT_PROCESS_TIME";
    /**
     * Name of CSV file with the results of the performance test of the count fingerprints.
     */
    private static final String CSV_COUNT_FINGERPRINT_PROCESS_RESULT = "CSV_COUNT_PROCESS_TIME";
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Initial capacity of the lists.
     */
    private final int initialCapacityValue = 200000;
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private  class variables">
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    /**
     * List which is contains all bit fingerprints
     */
    private ArrayList<ArrayList<String>> moleculeListforBitFingerprint = new ArrayList<>(this.initialCapacityValue);
    /**
     * Is a list that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private ArrayList<String> fragmentList = new ArrayList<>(this.initialCapacityValue);
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
     * Is a list in which all molecule fragments are stored that are read in from the CSV file.
     */
    private ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(this.initialCapacityValue);
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor.
     *
     * @param anArgs is the file name of the fragments (fragments are represented by unique SMILES in the file )
     *              that are read in to create the fingerprints.
     * @param anArgs2 is the file name of the molecules (molecules are represented by unique SMILES in the file )
     *               that are read in to create the fingerprints.
     * @param anArgs3 defines the size of the fragments to be read in as well as the molecule in the two files.
     * @throws IOException is thrown if the constructor is unable to open a text file for logging occurred exceptions.
     */
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
            this.resultsPrintWriter.println("Application initialized. Loading  files named " + anArgs + " and "+ anArgs2+ ".");
            this.resultsPrintWriter.println();
            // bit fingerprints process file
            File tmpBitFingerprintsResultFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_BIT_FINGERPRINT_PROCESS_RESULT + tmpProcessingTime + ".csv");
            FileWriter tmpBitFingerprintsResultWriter = new FileWriter(tmpBitFingerprintsResultFile, false);
            this.bitPrintWriter = new PrintWriter(tmpBitFingerprintsResultWriter);
            // count fingerprints process file
            File tmpCountFingerprintsResultFile = new File(this.workingPath + "/Results/" + PerformanceTest.CSV_COUNT_FINGERPRINT_PROCESS_RESULT + tmpProcessingTime + ".csv");
            FileWriter tmpCountFingerprintResultWriter = new FileWriter(tmpCountFingerprintsResultFile, false);
            this.countPrintWriter = new PrintWriter(tmpCountFingerprintResultWriter);
            // read in CSV files that contains fragments
            this.generateFingerprints(Integer.parseInt(anArgs3),1,2); // TODO delete arguments
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
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * The two text files (fragment and molecule file)
     * are read in and the data provided so that in the next step the
     * fingerprints can be calculated based on the information in these text files.
     *
     * @throws IOException is thrown if an error occurs when reading in the two text files.
     */
    private void LoadData() throws IOException {
        BufferedReader tmpFragmentSetReader;
        BufferedReader tmpMoleculeFragmentsReader;
        try {
            tmpFragmentSetReader =  new BufferedReader(new FileReader(this.moleculeFile));
            tmpMoleculeFragmentsReader = new BufferedReader(new FileReader(this.fragmentFile));
        } catch(IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("File is not readable");
        }
        //Read CSV file (fragmentation tab) to obtain fragments used to create fingerprints
        String tmpLine;
        String tmpSeparatorComma = ";";  // TODO separator as argument?
        this.fragmentList = new ArrayList<>();
        try {
            while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
                this.fragmentList.add(tmpSmilesOfFragments[0]);
            }
            this.fragmentList.remove(0);
        } catch(IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid fragment file. At least one line is not readable.");
        }
        // Read CSV file (itemization tab)
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        try {
            while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
                String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
                tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
            }
        } catch (IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid molecule file. At least one line is not readable");
        }
        tmpMoleculeFragmentsReader.close();
        List<String> tmpSeparateList;
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
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
                throw new IndexOutOfBoundsException("the line has not the right length");
            }
            this.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
            this.moleculeListforBitFingerprint.add(dataForGenerateBitFingerprint);
        }
    }
    //
    /**
     * Starts the generation of the fingerprints. And writes the results into the corresponding text files.
     *
     * @param aNumberOfMoleculesInProcess specifies how many fingerprints are to be generated.
     * @param aEndTime System time in ms.
     * @param aStartTime System time in ms.
     * @throws Exception is thrown if an error occurs when generating the fingerprints.
     */
    private void generateFingerprints(int aNumberOfMoleculesInProcess, long aEndTime, long aStartTime) throws Exception {
        try {
            this.LoadData();
        } catch (IOException anException) {
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
        this.resultsPrintWriter.println("\n\tGenerate bit fingerprints");
        this.resultsPrintWriter.println();
        this.bitPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListforBitFingerprint.size());
        this.bitPrintWriter.println("Number of molecules in process, Bit fingerprint process time in ms");
        for (int i = aNumberOfMoleculesInProcess; i <= this.moleculeListforBitFingerprint.size(); i+= aNumberOfMoleculesInProcess) {
            List<ArrayList<String>> tmpNumberOfMoleculesInProcess = this.moleculeListforBitFingerprint.subList(0, i);
            try {
             aStartTime = System.currentTimeMillis();
                for (ArrayList<String> tmpMolecule : tmpNumberOfMoleculesInProcess) {
                    IBitFingerprint bit = printer.getBitFingerprint(tmpMolecule);
                }
                 aEndTime = System.currentTimeMillis();
            } catch (Exception anException) {
                this.exceptionsPrintWriter.println("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list or molecule list.");
                this.appendToLogfile(anException);
                throw new Exception("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list or molecule list.");
            }
            resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcess.size() + " valid molecules.");
            resultsPrintWriter.println("Bit fingerprint generation took: " + (aEndTime - aStartTime) + " ms.");
            bitPrintWriter.println(tmpNumberOfMoleculesInProcess.size() + "," + (aEndTime - aStartTime));
        }
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("#########################################################################");
        this.resultsPrintWriter.println("\n\tGenerate count fingerprints");
        this.resultsPrintWriter.println();
        this.countPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.moleculeListforBitFingerprint.size());
        this.countPrintWriter.println("Number of molecules in process, Count fingerprint process time in ms.");
        for(int i = aNumberOfMoleculesInProcess; i<= this.moleculeFragmentList.size(); i+= aNumberOfMoleculesInProcess) {  //i+=2000)
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
    //</editor-fold>
}
