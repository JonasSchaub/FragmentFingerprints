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
 * @version 1.0.0.0
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
    private static final String CSV_BIT_FINGERPRINT_PROCESS_RESULT_FILE_NAME = "CSV_BIT_PROCESS_TIME";
    /**
     * Name of CSV file with the results of the performance test of the count fingerprints.
     */
    private static final String CSV_COUNT_FINGERPRINT_PROCESS_RESULT_FILE_NAME = "CSV_COUNT_PROCESS_TIME";
    /**
     * Name of the CSV file with the results of the generated bit fingerprints.
     */
    private static final String BIT_FINGERPRINT_RESULT_FILE_NAME = "BIT_FINGERPRINT_TIME";
    /**
     * Name of the CSV file with the results of the generated count fingerprints.
     */
    private static final String COUNT_FINGERPRINT_RESULT_FILE_NAME = "COUNT_FINGERPRINT_TIME";
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Initial capacity of the lists in which the data for generating the fingerprints is stored.
     */
    private final int INITIAL_CAPACITY_VALUE_NUMBER_OF_MOLECULES = 13000;
    /**
     * Initial capacity of the lists in which the data for generating the fingerprints is stored.
     */
    private final int INITIAL_CAPACITY_VALUE_NUMBER_OF_FRAGMENTS = 3988;
    /**
     * Value for determining the initial capacity of cards
     */
    private final int INITIAL_CAPACITY_VALUE = Math.round((4/3) + 1);
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
     * List contains all fragments of the molecules. These fragments are used in the creation of the bit fingerprint.
     */
    private ArrayList<ArrayList<String>> listOfMoleculeFragmentsList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_MOLECULES);
    /**
     * List that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private ArrayList<String> fragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_FRAGMENTS);
    /**
     * Stores the name/ID of the molecules from moleculeFile
     */
    private ArrayList<String> listOfMoleculeNames;
    /**
     * Input CSV file containing all fragments that and their frequencies for each given SMILES (molecule).
     */
    private File moleculeFile;
    /**
     * CSV file containing all fragments that are used for initialising the fingerprint.
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
    private PrintWriter bitFingerprintPrintWriter;
    /**
     * PrintWriter to log the generated count arrays.
     */
    private PrintWriter countFingerprintPrintWriter;
    /**
     * list stores each molecule as HashMap. The HashMap assigns
     * the unique SMILES (fragments) of the molecule to the frequency of the molecule,
     * which is read out from the molecule file.
     */
    private ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_MOLECULES);
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor.
     * All files listed below as parameters must be located in the same directory as the application JAR file.
     * In the constructor all necessary output files are created.
     * There are 3 different files to document the result of the performance snapshot. One file is created that
     * is essentially a console output and two others that contain the processing time required at each step in
     * a structured form (CSV). Another file contains the generated fingerprints for each molecule.
     * The constructor executes the complete application. Everything the application is supposed to do is
     * realized in the constructor, including the shutdown at the end of the application.
     *
     * @param anArgs the command line arguments, anArgs[0] and anArgs[1] must be the names of the text files to load.
     *               The text files must be located in the same directory as the application's JAR file. The command
     *               line anArgs[3] must be the name of the path, where the performance test results are
     *               stored. The command line anArgs[2] must be the number of the molecules in process.
     *               <ul>
     *                <li>anArgs[0]:
     *                              <ul>
     *                          <li> is a semicolon separated CSV file with header</li>
     *                          <li>First line must contain unique SMILES strings.
     *                          These are set as predefined fragments when the fingerprint is initialized.</li>
     *                          <li>All other lines in the file also the header are ignored.</li>
     *                      </ul>
     *                </li>
     *                <li>anArgs[1]:
     *                        <ul>
     *                           <li> is a semicolon-separated CSV file with header line</li>
     *                           <li> Each row shows a molecule with its fragments</li>
     *                           <li> first column: Molecule name/ID</li>
     *                           <li>second column: unique SMILES of the molecule</li>
     *                           <li>third and following rows: alternating pairs of fragment SMILES and abundances
     *                               in the respective molecule</li>
     *                           <li>In this file, the header and columns 1 and 2 are ignored.</li>
     *                           <li>The SMILES fragments from this file are matched with the predefined fragments
     *                               to create the fingerprint.</li>
     *                </ul>
     *               </li>
     *               </ul>
     *
     *
     * @throws IOException is thrown if the constructor is unable to open a text file for logging occurred exceptions.
     * @throws IllegalArgumentException is thrown if the given arguments are invalid.
     */
    public PerformanceTest(String[] anArgs) throws IOException, IllegalArgumentException {
        if(anArgs.length !=4) {
            throw new IllegalArgumentException("Four arguments (three file names and the number of the molecules in process) are required.");
        }
        String tmpInputFilePath =  (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        /*Set up exception log file*/
        if(anArgs[2].isEmpty() || anArgs[2].isBlank()) {
            throw new IllegalArgumentException("The 3rd pathname is empty or blank.");
        }
        this.workingPath = (new File(anArgs[2]).getAbsoluteFile().getAbsolutePath()) + File.separator;
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpOutputPath = this.workingPath + "Results" + File.separator;
        String tmpProcessingTime = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        new File(tmpOutputPath).mkdirs();
        File tmpExceptionsLogFile = new File( tmpOutputPath
                + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME + tmpProcessingTime);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        this.exceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        this.exceptionsPrintWriter.println("#########################################################################");
        this.exceptionsPrintWriter.println("Processing Time: " + tmpProcessingTime);
        this.exceptionsPrintWriter.println();
        this.exceptionsPrintWriter.flush();
        try {
            // Load fragment files
            this.moleculeFile = new File(tmpInputFilePath + anArgs[0]);
            this.fragmentFile = new File(tmpInputFilePath + anArgs[1]);
            // test the input files for validity
            try (FileInputStream tmpMoleculeFileInputStream = new FileInputStream(this.moleculeFile);
                 FileInputStream tmpFragmentFileInputStream = new FileInputStream(this.fragmentFile)) {
            } catch (FileNotFoundException | SecurityException anException) {
                this.appendToLogfile(anException);
                throw new IllegalArgumentException("One or more files (name) are invalid: " + anException.getMessage());
            }
            // results files
            File tmpResultsLogFile = new File(tmpOutputPath + PerformanceTest.RESULTS_FILE_NAME + tmpProcessingTime + ".txt");
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            this.resultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            this.resultsPrintWriter.println("#########################################################################");
            this.resultsPrintWriter.println();
            this.resultsPrintWriter.println("Time: " + tmpProcessingTime);
            this.resultsPrintWriter.println();
            this.resultsPrintWriter.println("Application initialized. Loading  files named " + anArgs[0] + " and "+ anArgs[1] + ".");
            this.resultsPrintWriter.println();
            // bit fingerprints process file
            File tmpBitFingerprintsOutputFile = new File(tmpOutputPath + PerformanceTest.CSV_BIT_FINGERPRINT_PROCESS_RESULT_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpBitFingerprintsResultWriter = new FileWriter(tmpBitFingerprintsOutputFile, false);
            this.bitPrintWriter = new PrintWriter(tmpBitFingerprintsResultWriter);
            // count fingerprints process file
            File tmpCountFingerprintsOutputFile = new File(tmpOutputPath + PerformanceTest.CSV_COUNT_FINGERPRINT_PROCESS_RESULT_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpCountFingerprintResultWriter = new FileWriter(tmpCountFingerprintsOutputFile, false);
            this.countPrintWriter = new PrintWriter(tmpCountFingerprintResultWriter);
            // bit fingerprint result file
            File tmpBitArrayFingerprintResultFile = new File(tmpOutputPath + PerformanceTest.BIT_FINGERPRINT_RESULT_FILE_NAME + tmpProcessingTime+ ".txt");
            FileWriter tmpBitArrayFingerprintResultWriter = new FileWriter(tmpBitArrayFingerprintResultFile, false);
            this.bitFingerprintPrintWriter = new PrintWriter(tmpBitArrayFingerprintResultWriter);
            // count fingerprint result file
            File tmpCountArrayFingerprintResultFile = new File(tmpOutputPath + PerformanceTest.COUNT_FINGERPRINT_RESULT_FILE_NAME + tmpProcessingTime+ ".txt");
            FileWriter tmpCountArrayFingerprintResultWriter = new FileWriter(tmpCountArrayFingerprintResultFile, false);
            this.countFingerprintPrintWriter = new PrintWriter(tmpCountArrayFingerprintResultWriter);
            // read in CSV files that contain fragments
            try {
                this.importDataFromTextFile();
            } catch (IOException anException) {
                this.exceptionsPrintWriter.println("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in.");
                this.exceptionsPrintWriter.flush();
                this.appendToLogfile(anException);
                throw new Exception("Fragment load ERROR. Unsuitable (structure) CSV files were tried to be read in. ");
            }
            // generate bit and count fingerprints
            this.generateFingerprints(Integer.parseInt(anArgs[3]));
            this.resultsPrintWriter.flush();
            this.bitPrintWriter.flush();
            this.countPrintWriter.flush();
            this.bitFingerprintPrintWriter.flush();
            this.countFingerprintPrintWriter.flush();
            this.exceptionsPrintWriter.flush();
            System.out.println("Application is finished");
        } catch (Exception anException) {
            this.appendToLogfile(anException);
            anException.printStackTrace(System.err);
            System.exit(1);
        } finally {
            this.exceptionsPrintWriter.close();
            this.resultsPrintWriter.close();
            this.bitFingerprintPrintWriter.close();
            this.countPrintWriter.close();
            this.countFingerprintPrintWriter.close();
            this.bitPrintWriter.close();
        }
    }
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * The two text files (fragment and molecule file)
     * are read in and the data provided so that in the next step the
     * fingerprints can be calculated based on the information in these text files.
     * For the method to work properly, it is important that the two files are in a specific form defined
     * for the method. Thus, each line in the fragment file corresponds to a fragment. Components
     * within a line must be separated by
     * a semicolon. The fragment represented by a unique SMILES must be at the very front (first position) of the line.
     * All other information within the line is ignored.
     * The same applies to the molecule file. Here, too, a molecule is defined in each line. However, each line
     * starts with the name/ID of the molecule and the SMILES encoding corresponding to the molecule structure.
     * Then, the fragments of the molecule and the corresponding frequency of the fragment in the molecule follow
     * alternately in the line. Again, all components are separated by a semicolon.
     *
     * @throws IOException is thrown if an error occurs when reading in the two text files.
     */
    private void importDataFromTextFile() throws IOException {
        try (BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader(this.moleculeFile));
             BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader(this.fragmentFile))) {
            //Read CSV file (fragments file) to obtain fragments used to create fingerprints
            String tmpLine;
            this.fragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_FRAGMENTS);
            try {
                while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                    String[] tmpSmilesOfFragments = tmpLine.split(","); // this.LINE_SEPARATOR_SEMICOLON
                    this.fragmentList.add(tmpSmilesOfFragments[0]);
                }
                // removing header line value
                this.fragmentList.remove(0);
            } catch(IOException anException) {
                this.appendToLogfile(anException);
                throw new IOException("invalid fragment file. At least one line is not readable.");
            }
            // Read CSV file (molecules file)
            List<List<String>> tmpListOfMoleculesFragmentsAndFrequenciesList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_MOLECULES);
            String tmpMoleculeLine;
            try {
                while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
                    String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(this.LINE_SEPARATOR_SEMICOLON);
                    tmpListOfMoleculesFragmentsAndFrequenciesList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
                }
            } catch (IOException anException) {
                this.appendToLogfile(anException);
                throw new IOException("invalid molecule file. At least one line is not readable");
            }
            List<String> tmpMoleculeFragmentsAndFrequenciesList;
            this.listOfMoleculeNames = new ArrayList<>();
            for (int tmpCurrentLineIndex = 1; tmpCurrentLineIndex < tmpListOfMoleculesFragmentsAndFrequenciesList.size(); tmpCurrentLineIndex++) {
                tmpMoleculeFragmentsAndFrequenciesList = tmpListOfMoleculesFragmentsAndFrequenciesList.get(tmpCurrentLineIndex);
                this.listOfMoleculeNames.add(tmpMoleculeFragmentsAndFrequenciesList.get(0));
                List<String> tmpListWithoutNameAndMoleculeSmiles = tmpMoleculeFragmentsAndFrequenciesList.subList(2, tmpMoleculeFragmentsAndFrequenciesList.size());
                HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>(tmpListWithoutNameAndMoleculeSmiles.size()*INITIAL_CAPACITY_VALUE);
                ArrayList<String> tmpDataToGenerateBitFingerprint = new ArrayList<>();
                try {
                    for (int i = 0; i < tmpListWithoutNameAndMoleculeSmiles.size(); i++) {
                        if (i % 2 == 0) { // magic number to store the fragment SMILES and their frequencies from the file into a HashMap
                            tmpMoleculeFragmentsMap.put(tmpListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(tmpListWithoutNameAndMoleculeSmiles.get(i + 1)));
                            tmpDataToGenerateBitFingerprint.add(tmpListWithoutNameAndMoleculeSmiles.get(i));
                        }
                    }
                } catch (IndexOutOfBoundsException anException) {
                    this.appendToLogfile(anException);
                    throw new IndexOutOfBoundsException("The line has not the right length");
                }
                this.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
                this.listOfMoleculeFragmentsList.add(tmpDataToGenerateBitFingerprint);
            }
        } catch (IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("File is not readable");
        }
    }

    /*
        //Read CSV file (fragments file) to obtain fragments used to create fingerprints
        String tmpLine;
        this.fragmentList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_FRAGMENTS);
        try {
            while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                String[] tmpSmilesOfFragments = tmpLine.split(this.LINE_SEPARATOR_SEMICOLON);
                this.fragmentList.add(tmpSmilesOfFragments[0]);
            }
            // removing header line value
            this.fragmentList.remove(0);
        } catch(IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid fragment file. At least one line is not readable.");
        }
        finally {
            tmpFragmentSetReader.close();
        }
        // Read CSV file (molecules file)
        List<List<String>> tmpListOfMoleculesFragmentsAndFrequenciesList = new ArrayList<>(this.INITIAL_CAPACITY_VALUE_NUMBER_OF_MOLECULES);
        String tmpMoleculeLine;
        try {
            while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
                String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(this.LINE_SEPARATOR_SEMICOLON);
                tmpListOfMoleculesFragmentsAndFrequenciesList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
            }
        } catch (IOException anException) {
            this.appendToLogfile(anException);
            throw new IOException("invalid molecule file. At least one line is not readable");
        }
        finally {
            tmpMoleculeFragmentsReader.close();
        }
        List<String> tmpMoleculeFragmentsAndFrequenciesList;
        this.listOfMoleculeNames = new ArrayList<>();
        for (int tmpCurrentLineIndex = 1; tmpCurrentLineIndex < tmpListOfMoleculesFragmentsAndFrequenciesList.size(); tmpCurrentLineIndex++) {
            tmpMoleculeFragmentsAndFrequenciesList = tmpListOfMoleculesFragmentsAndFrequenciesList.get(tmpCurrentLineIndex);
            this.listOfMoleculeNames.add(tmpMoleculeFragmentsAndFrequenciesList.get(0));
            List<String> ListWithoutNameAndMoleculeSmiles = tmpMoleculeFragmentsAndFrequenciesList.subList(2, tmpMoleculeFragmentsAndFrequenciesList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            ArrayList<String> tmpDataToGenerateBitFingerprint = new ArrayList<>();
            try {
                for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                    if (i % 2 == 0) { // magic number to store the data from the file into a HashMap
                        tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                        tmpDataToGenerateBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                    }
                }
            } catch (IndexOutOfBoundsException anException) {
                this.appendToLogfile(anException);
                throw new IndexOutOfBoundsException("The line has not the right length");
            }
            this.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
            this.listOfMoleculeFragmentsList.add(tmpDataToGenerateBitFingerprint);
        }

     */
    //
    /**
     * Starts the generation of the fingerprints. And writes the results into the corresponding text file.
     *
     * @param aNumberOfMoleculesInProcess bin size that specifies how many fingerprints to generate in one iteration.
     * @throws Exception is thrown if an error occurs when generating the fingerprints.
     */
    private void generateFingerprints(int aNumberOfMoleculesInProcess) throws Exception {
        FragmentFingerprinter tmpFragmentFingerprinter = new FragmentFingerprinter(this.fragmentList);
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("Number of molecules: " + this.listOfMoleculeFragmentsList.size());
        System.out.println("Number of molecules: " + this.listOfMoleculeFragmentsList.size());
        this.resultsPrintWriter.println("Number of fragments: " + this.fragmentList.size());
        System.out.println("Number of fragments: " + this.fragmentList.size());
        this.resultsPrintWriter.println("Completed loading of molecules and fragments");
        System.out.println("Completed loading of molecules and fragments");
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("\n\tGenerate bit fingerprints");
        this.resultsPrintWriter.println();
        this.bitPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.listOfMoleculeFragmentsList.size());
        this.bitPrintWriter.println("Molecules processed, Bit fingerprint process time in ms");
        this.bitFingerprintPrintWriter.println("Molecule name/ID, bit fingerprint");
        this.countFingerprintPrintWriter.println("Molecule name/ID, count fingerprint");
        for (int i = aNumberOfMoleculesInProcess; i <= this.listOfMoleculeFragmentsList.size(); i+= aNumberOfMoleculesInProcess) {
            List<ArrayList<String>> tmpNumberOfMoleculesInProcess = this.listOfMoleculeFragmentsList.subList(0, i);
            long tmpStartTime = System.currentTimeMillis();
            for (ArrayList<String> tmpMolecule : tmpNumberOfMoleculesInProcess) {
                try {
                    IBitFingerprint tmpBitFingerprint = tmpFragmentFingerprinter.getBitFingerprint(tmpMolecule);
                } catch (Exception anException) {
                    this.exceptionsPrintWriter.println("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                    this.appendToLogfile(anException);
                    throw new Exception("Bit fingerprint generation ERROR. There may be incorrect/invalid elements in the fragment list oder molecule list");
                }
            }
            long tmpEndTime = System.currentTimeMillis();
            this.resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcess.size() + " valid molecules.");
            this.resultsPrintWriter.println("Bit fingerprint generation took: " + (tmpEndTime - tmpStartTime) + " ms.");
            this.bitPrintWriter.println(tmpNumberOfMoleculesInProcess.size() + "," + (tmpEndTime - tmpStartTime));
        }
        int[] tmpBitArray;
        for(Iterator tmpMoleculeNameIterator = this.listOfMoleculeNames.iterator(), tmpMolecule = this.listOfMoleculeFragmentsList.iterator(); tmpMoleculeNameIterator.hasNext() && tmpMolecule.hasNext();) {
            List<String> tmpListOfMolecule = (List<String>) tmpMolecule.next();
            String tmpMoleculeNameOrID = (String) tmpMoleculeNameIterator.next();
            tmpBitArray = tmpFragmentFingerprinter.getBitArray(tmpListOfMolecule);
            this.bitFingerprintPrintWriter.println(tmpMoleculeNameOrID + "," + java.util.Arrays.toString(tmpBitArray));
        }
        this.resultsPrintWriter.println();
        this.resultsPrintWriter.println("#########################################################################");
        this.resultsPrintWriter.println("\n\tGenerate count fingerprints");
        this.resultsPrintWriter.println();
        this.countPrintWriter.println("Number of fragments: " + this.fragmentList.size() + " and number of molecules: " + this.listOfMoleculeFragmentsList.size());
        this.countPrintWriter.println("Number of molecules in process, Count fingerprint process time in ms.");
        for(int i = aNumberOfMoleculesInProcess; i<= this.moleculeFragmentList.size(); i+= aNumberOfMoleculesInProcess) {
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
            long tmpEndTime = System.currentTimeMillis();
            this.resultsPrintWriter.println("Processing " + tmpNumberOfMoleculesInProcessList.size() + " valid molecules.");
            this.resultsPrintWriter.println("Count fingerprint generation took: " + (tmpEndTime - start) + " ms.");
            this.resultsPrintWriter.println();
            this.countPrintWriter.println(tmpNumberOfMoleculesInProcessList.size() + "," + (tmpEndTime - start));
        }
        System.out.println(this.moleculeFragmentList + "-----molecule fragment list");

        int[] tmpCountArray;
        /*
        for(Iterator tmpMoleculeNameIterator = this.listOfMoleculeNames.iterator(), tmpMolecule = this.moleculeFragmentList.iterator(); tmpMoleculeNameIterator.hasNext() && tmpMolecule.hasNext();) {
            HashMap<String, Integer> tmpListOfMolecule = (HashMap<String, Integer>) tmpMolecule.next();
            System.out.println(tmpListOfMolecule+ "-----list of molecule");
            String tmpMoleculeNameOrID = (String) tmpMoleculeNameIterator.next();
            tmpCountArray = tmpFragmentFingerprinter.getCountArray(tmpListOfMolecule);
            System.out.println(java.util.Arrays.toString(tmpCountArray)+ "---tmpCountARray");
            this.countFingerprintPrintWriter.println(tmpMoleculeNameOrID + "," + java.util.Arrays.toString(tmpCountArray));


        }

         */
        for(HashMap<String, Integer> map : this.moleculeFragmentList) {
            System.out.println(map+ "---map");
            System.out.println(java.util.Arrays.toString(tmpFragmentFingerprinter.getCountArray(map)) + "-----count array");
        }

        System.out.println("Bit and count fingerprints were generated successfully.");
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
        StringWriter tmpStringWriter = new StringWriter();
        anException.printStackTrace(new PrintWriter(tmpStringWriter));
        String tmpStackTrace = tmpStringWriter.toString();
        this.exceptionsPrintWriter.println(tmpStackTrace);
        if (this.exceptionsPrintWriter != null) {
            this.exceptionsPrintWriter.println();
            this.exceptionsPrintWriter.flush();
        }
    }
    //</editor-fold>
}
