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

package de.unijena.cheminf.fragment.fingerprint;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Class to generate fragment fingerprints. Bit and count fragment fingerprints can be generated.
 * Fragment fingerprints are key-based fingerprints.
 * Thus, the class requires predefined structures/fragments in
 * the form of unique SMILES to create the fingerprint. These structures must be passed when the
 * class is instantiated (in the constructor). The class implements the interface IFragmentFingerprinter,
 * which inherits the IFingerprinter (CDK), which allows the class to compute fingerprints in 2 ways.
 * The first way to calculate a bit or count fingerprint is to perform a substructure comparison with all
 * predefined fragments for a given IAtomContainer. The fingerprint created by the substructure search is based on
 * the CDK class SubstructureFingerprinter. The predefined fragment SMILES are interpreted as SMARTS patterns by the
 * SubstructureFingerprinter class. The second way to calculate fingerprints is by comparing
 * given fragments, which are in the form of unique SMILES, with the predefined fragments.
 * The second possibility is thus based on a pure comparison of strings. It is important to note that the two
 * different ways of creating fingerprints can produce different results.
 *
 * @author Betuel Sevindik
 * @version 1.0.0.0
 */
public class FragmentFingerprinter implements IFragmentFingerprinter {
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * The array containing all the unique predefined (key) SMILES fragments based
     * on which the fingerprints are then created. This set of fragment/unique SMILES
     * must be included when initializing this class.
     */
    private final String[] fragmentFingerprintMasterVectorArray;
    /**
     * The fragmentArray is converted into a HashMap to speed up the matching of the unique SMILES.
     * The Map maps the unique SMILES of the predefined fragments to the position they have in the array.
     */
    private final HashMap<String, Integer> uniqueSmilesToPositionMap;
    //</editor-fold>
    //
    //<editor-fold desc="private static final class variables" defaultstate="collapsed">
    /**
     * Version of fragment fingerprinter
     */
    private static final String FRAGMENT_FINGERPRINTER_VERSION = "1.0.0.0";
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    /**
     * Bit set fingerprint for storing the calculated fragment bit fingerprint.
     */
    private BitSetFingerprint cacheBitSetFingerprint;
    /**
     * Count fingerprint for storing the calculated fragment count fingerprint.
     */
    private CountFingerprint cacheCountFingerprint;
    //setting booleans
    /**
     * Private boolean for setting 'Use Cache Fingerprint'. This setting determines whether an already created, cached
     * fingerprint should be used in internal calculations when a second fingerprint calculation is started.
     */
    private boolean cacheFingerprintSetting = false;
    //</editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor.
     * Initialization of the fragment fingerprinter by using a user-defined
     * set of fragments in the form of unique SMILES.
     * If the list passed during initialization contains duplicates, they will be removed.
     * The number of predefined fragments specified by the user may then differ from the actual number of
     * key fragments present, as duplicates are removed. This means that duplicate fragment SMILES strings in the input
     * list are ignored and are not part of the fingerprint multiple times.
     *
     * @param aFragmentsForMasterVectorList is the ist in which the predefined fragments are stored.
     * @throws NullPointerException is thrown if the list aFragmentsForMasterVectorList is null.
     * @throws IllegalArgumentException is thrown if the list contains blank strings.
     */
    public FragmentFingerprinter(List<String> aFragmentsForMasterVectorList) throws NullPointerException, IllegalArgumentException {
        // Check whether aFragmentsForMasterVectorList is null or whether there are elements (strings) in the list that are empty.
        this.validityCheckOfParameterList(aFragmentsForMasterVectorList,"aFragmentsForMasterVectorList (list of string instances) is null.",
                "aFragmentsForMasterVectorList (at least one list element) is null.",
                "aFragmentsForMasterVectorList (at least one list element) is blank/empty.");
        this.fragmentFingerprintMasterVectorArray = aFragmentsForMasterVectorList.toArray(new String[aFragmentsForMasterVectorList.size()]);
        this.uniqueSmilesToPositionMap = this.buildUniqueSmilesToPositionMap(this.fragmentFingerprintMasterVectorArray);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Overriden public methods">
    /**
     * Method to generate the bit fingerprint.
     * An entered list of unique SMILES is compared with the predefined fragments.
     * If there is a match, the position of the unique SMILES is determined from the map and set to true in the
     * initialized BitSet. The method is intended to generate the fingerprint for one molecule but can in principle be
     * applied to any fragment set, e.g. originating from a cluster of multiple molecules.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * To be able to calculate the fingerprint for a molecule, the fragments should belong to one molecule.
     * @return BitSet. BitSet is a CDK class that implements the IBitFingerprint interface of CDK.
     * This allows methods to be used that return useful information from the calculated bit fingerprint,
     * such as the number of positive bits in the fingerprint, etc.
     * @throws NullPointerException is thrown if the list aListOfUniqueSmiles is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    @Override
    public IBitFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"Given list of string instances is null.",
                "Given list includes at least one null element.",
                "Given list includes at least one blank/empty element.");
        if (this.cacheFingerprintSetting && this.cacheBitSetFingerprint != null) {
            return this.cacheBitSetFingerprint;
        }
        BitSet tmpBitSet = new BitSet(this.uniqueSmilesToPositionMap.size());
        if (this.cacheFingerprintSetting && this.cacheCountFingerprint != null) {
            for (int i = 0; i < this.cacheCountFingerprint.size(); i++) {
                if (this.cacheCountFingerprint.getCount(i) > 0) {
                    tmpBitSet.set(i, true);
                }
            }
        } else {
            for (String tmpSmiles : aListOfUniqueSmiles) {
                if (this.uniqueSmilesToPositionMap.containsKey(tmpSmiles)) {
                    tmpBitSet.set(this.uniqueSmilesToPositionMap.get(tmpSmiles), true);
                }
            }
        }
        if (this.cacheFingerprintSetting) {
            this.cacheBitSetFingerprint = new BitSetFingerprint(tmpBitSet);
            return this.cacheBitSetFingerprint;
        } else {
            return new BitSetFingerprint(tmpBitSet);
        }
    }
    //
    //<editor-fold desc="Public Methods - float variants">
    /**
     * Public method for creating a float "bit set fingerprint" in form of a float[] array.
     * For this, a given list of SMILES is compared with the pre-defined fragments (master vector of fragments for fingerprint definition).
     * On the basis of the given SMILES to position map, each SMILES is checked whether the defined fingerprint includes the fragment.
     * If the defined fingerprint contains the fragment, it's respective position inside the array is set to '1.0f'.
     *
     * @param aFragmentsUniqueSmilesList to be checked against the pre-defined fingerprint
     * @param aPreInitFloatArray         pre-initialized float[] array with size of pre-defined fragment fingerprint and additional space for descriptive float components
     */
    public void getFloatBitFingerprint(
            //aFragmentsUniqueSmilesList contains fragments of ONE molecule
            List<String> aFragmentsUniqueSmilesList,
            float[] aPreInitFloatArray) {
        for (String tmpSmiles : aFragmentsUniqueSmilesList) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpSmiles)) {
                aPreInitFloatArray[this.uniqueSmilesToPositionMap.get(tmpSmiles)] = 1.0f;
            }
        }
    }
    //
    /**
     * Public method for generating a float count fingerprint in form of a float[] array.
     * This fingerprint projects the frequency (-> count) of occurrence of the pre-defined fingerprint fragments within
     * the given SMILES.
     * <p>
     * For example: The pre-defined fragment fingerprint contains the SMILES for methane as "C", for an ether substructure
     * as "*O*", and for a hydroxide group as "[H]OC". The given SMILES list contains only three "C" fragments. Then,
     * the generated fingerprint would look like the following: [3.0f, 0.0f, 0.0f]
     * </p>
     *
     * @param aFragmentsUniqueSmilesList of fragment SMILES to be compared to the pre-defined fingerprint
     * @param aPreInitFloatArray         pre-initialized float[] in which the fingerprint is to be stored
     */
    public void getFloatCountFingerprint(
            List<String> aFragmentsUniqueSmilesList,
            float[] aPreInitFloatArray) {
        for (String tmpSmiles : aFragmentsUniqueSmilesList) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpSmiles)) {
                aPreInitFloatArray[this.uniqueSmilesToPositionMap.get(tmpSmiles)]++;
            }
        }
    }
    //
    /**
     * Public method to get a float[][] matrix containing the fingerprints for the given SMILES in regards to the
     * pre-defined fragment fingerprint.
     * It shall be noted that each List in the given List (List of Lists) represents ONE molecule's fragments.
     * Therefor, the whole list represents the entirety of fragments.
     * Further, a setting whether to use "bit set" or "count/frequency of fragment" is available.
     *
     * @param aFragmentsUniqueSmilesListsArrayList with Lists of SMILES to compare with the defined fingerprint
     * @param aFloatDataMatrix to be filled with the float fingerprints of the given fragments
     * @param anUseBitArrayStatement setting whether "bit set" or "count/frequency" should be used for the matrix
     */
    public void getFragmentsComponentsFloatMatrix(
            List<List<String>> aFragmentsUniqueSmilesListsArrayList,
            float[][] aFloatDataMatrix,
            boolean anUseBitArrayStatement
    ) {
        if (aFragmentsUniqueSmilesListsArrayList.size() <= aFloatDataMatrix.length) {
            //float[row count][column count]
            for (int i = 0; i < aFloatDataMatrix.length; i++) {
                if (anUseBitArrayStatement) {
                    this.getFloatBitFingerprint(aFragmentsUniqueSmilesListsArrayList.get(i), aFloatDataMatrix[i]);
                } else {
                    this.getFloatCountFingerprint(aFragmentsUniqueSmilesListsArrayList.get(i), aFloatDataMatrix[i]);
                }
            }
        }
    }
    //</editor-fold>
    //
    /**
     * Generates count fingerprint for a molecule based on its fragments represented by unique SMILES strings
     * in the key set and their frequencies in the value set of the given map. Given fragment SMILES codes that
     * are not part of the set given at initialisation of this class, are ignored. The frequencies of those matching
     * with the predefined set are used to construct the fingerprint. The method is intended to generate the fingerprint
     * for one molecule but can in principle be applied to any fragment set, e.g. originating from a cluster of
     * multiple molecules
     *
     * @param aSmilesToFrequencyMap map usually represents a molecule by representing the fragments of
     * the molecule by unique SMILES in the key set and indicating their frequency in the value set. In principle,
     * however,such a map can be applied to any set of fragments.
     * To be able to calculate the fingerprint for a molecule, the fragments must belong to a molecule.
     * @return count fingerprint
     * @throws NullPointerException  is thrown if the map aSmilesToFrequencyMap is
     * null or contains keys or values that are null respectively.
     * @throws IllegalArgumentException is thrown if the map aSmilesToFrequencyMap
     * contains keys or values that are blank/empty, respectively.
     */
    @Override
    public ICountFingerprint getCountFingerprint(Map<String, Integer> aSmilesToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        HashMap<Integer, Integer> tmpPositionToFrequencyMap = new HashMap<>(
                (int) (this.uniqueSmilesToPositionMap.size() * 1.5f),
                0.75f);
        Objects.requireNonNull(aSmilesToFrequencyMap, "Given map of string and integer instances is null.");
        for (Map.Entry<String, Integer> tmpEntry : aSmilesToFrequencyMap.entrySet()) {
            if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                throw new NullPointerException("Given map of string and integer instances contains " +
                        "instances that are null.");
            } else if (tmpEntry.getKey().isEmpty() || tmpEntry.getKey().isBlank()) {
                throw new IllegalArgumentException("Given map of strings an integer instances contains strings that are blank/empty.");
            } else if (this.uniqueSmilesToPositionMap.containsKey(tmpEntry.getKey())) {
                tmpPositionToFrequencyMap.put(this.uniqueSmilesToPositionMap.get(tmpEntry.getKey()), tmpEntry.getValue());
            }
        }
        return new CountFingerprint(this.fragmentFingerprintMasterVectorArray, tmpPositionToFrequencyMap);
    }
    //
    /**
     * Generates a count fingerprint for a molecule based on its fragments, represented by unique SMILES in the
     * list given as parameters. Given fragment SMILES codes that are not part of the set given at initialisation of
     * this class, are ignored. The frequencies of those matching with the predefined set are used to construct the
     * fingerprint. The frequency of individual fragments depends on how often they occur in the specified list.
     * Duplicates are thus allowed in this list.
     * The method is intended to generate the fingerprint for one molecule but can in principle be applied to any
     * fragment set, e.g. originating from a cluster of multiple molecules.
     *
     *
     * @param aUniqueSmilesList is a list that stores fragments in the form of unique SMILES.
     * If a fragment occurs more than once in the molecule, it is also present more than
     * once in the list. To be able to calculate the fingerprint for a molecule,
     * the fragments should belong to one molecule.
     * @return count fingerprint
     * @throws NullPointerException is thrown if the list aUniqueSmilesToFrequencyList is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    @Override
    public ICountFingerprint getCountFingerprint(List<String> aUniqueSmilesList) throws NullPointerException, IllegalArgumentException {
        HashMap<String, Integer> tmpUniqueSmilesToFrequencyCountMap = new HashMap<>((int) (this.uniqueSmilesToPositionMap.size() * 1.5f), 0.75f);
        Objects.requireNonNull(aUniqueSmilesList, "aUniqueSmilesToFrequencyList (list of string instances) is null.");
        for (String tmpSmiles : aUniqueSmilesList) {
            Objects.requireNonNull(tmpSmiles, "aUniqueSmilesToFrequencyList (at least one list element) is null.");
            if(tmpSmiles.isBlank()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyList (at least one list element) is blank/empty.");
            }
            if (!tmpUniqueSmilesToFrequencyCountMap.containsKey(tmpSmiles)) {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, 1);
            } else {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, tmpUniqueSmilesToFrequencyCountMap.get(tmpSmiles) + 1);
            }
        }
        return this.getCountFingerprint(tmpUniqueSmilesToFrequencyCountMap);
    }
    //
    /**
     * {@inheritDoc}
     */
    @Override
    public String getVersionDescription() {
        return getClass().getSimpleName() + "/" + FragmentFingerprinter.FRAGMENT_FINGERPRINTER_VERSION +
                ' ' + "num_bits" + '=' + this.uniqueSmilesToPositionMap.size();
    }
    //
    /**
     * {@inheritDoc}
     */
    @Override
    public BitSet getFingerprint(IAtomContainer mol) throws CDKException {
       IBitFingerprint tmpAtomContainerBasedBitFingerprint =  this.getBitFingerprint(mol);
       return tmpAtomContainerBasedBitFingerprint.asBitSet();
    }
    //
    /**
     * {@inheritDoc}
     * @see SubstructureFingerprinter
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        String[] tmpPredefinedFragments = this.getPredefinedFragmentArrayWithoutDuplicates();
        SubstructureFingerprinter tmpSubstructureFingerprint = new SubstructureFingerprinter(tmpPredefinedFragments);
        IBitFingerprint tmpBitFingerprintBySubstructureSearch = tmpSubstructureFingerprint.getBitFingerprint(container);
        return tmpBitFingerprintBySubstructureSearch;
    }
    //
    /**
     * {@inheritDoc}
     * @see SubstructureFingerprinter
     */
    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        String[] tmpPredefinedFragments = this.getPredefinedFragmentArrayWithoutDuplicates();
        SubstructureFingerprinter tmpSubstructureFingerprint = new SubstructureFingerprinter(tmpPredefinedFragments);
        ICountFingerprint tmpCountFingerprintBySubstructureSearch = tmpSubstructureFingerprint.getCountFingerprint(container);
        return tmpCountFingerprintBySubstructureSearch;
    }
    //
    /**
     * UnsupportedOperationException. This method is not supported.
     * {@inheritDoc}
     *
     * @throws UnsupportedOperationException method is not supported
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {
       throw new UnsupportedOperationException();
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since the FragmentFingerprinter is a key-based fingerprint, the size of the fingerprint is equal
     * to the number of predefined fragments (unique SMILES) if the list of key fragments passed during
     * initialization does not contain duplicates, otherwise the size of the fingerprint may be smaller
     * than the number of fragments passed since duplicates are removed.
     * Which means that duplicate fragment SMILES strings are ignored during initialization and are not
     * part of the fingerprint multiple times.
     *
     * @return int
     */
    @Override
    public int getSize() {
        return this.uniqueSmilesToPositionMap.size();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    //
    /**
     * Set method for 'Use Cache Fingerprint' setting. The setting determines if already created fingerprints should be used
     * in internal methods in a following calculation.
     *
     * @param aSettingValue boolean setting value
     */
    public void setCacheFingerprintSetting(boolean aSettingValue) {
        this.cacheFingerprintSetting = aSettingValue;
    }
    //
    /**
     * Returns the bit definitions i.e. which  bit stands for which fragment SMILES.
     * Important, the number of possible bit definitions may differ from the number of key
     * fragments passed during initialization, since duplicates are removed.
     *
     * @param aBit position in the fingerprint.
     * @return unique SMILES corresponding to the specified position.
     * @throws IllegalArgumentException is thrown if the given bit position is not present in the fingerprint.
     */
    public String getBitDefinition(int aBit) throws IllegalArgumentException {
        String[] tmpPredefinedFragmentsInArray = this.getPredefinedFragmentArrayWithoutDuplicates();
        if(aBit < tmpPredefinedFragmentsInArray.length && aBit >= 0) {
            return tmpPredefinedFragmentsInArray[aBit];
        } else {
            throw new IllegalArgumentException("This bit is not defined/present in the fingerprint.");
        }
    }
    //
    /**
     * Returns bit array for specified list.
     * The size of the array corresponds to the number of predefined (key) fragments passed during initialization.
     * However, the size may differ if there are duplicates in the specified predefined fragments, as they
     * will be ignored/removed.
     * This method is only available for bit fingerprints based on unique SMILES comparisons.
     *
     * @param aListOfUniqueSmiles is a list that stores molecule fragments or arbitrary fragments
     * in the form of unique SMILES.
     * @return int[] bit array
     * @throws NullPointerException is thrown if the list aListOfUniqueSmiles is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    public int[] getBitArray(List<String> aListOfUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        //converts BitSetFingerprint to BitArray
        int[] tmpReturnArray = new int[this.fragmentFingerprintMasterVectorArray.length];
        if (this.cacheFingerprintSetting && this.cacheBitSetFingerprint != null) {
            for (int tmpPositivePosition : this.cacheBitSetFingerprint.getSetbits()) {
                tmpReturnArray[tmpPositivePosition] = 1;
            }
        } else {
            for (int tmpPositivePosition : this.getBitFingerprint(aListOfUniqueSmiles).getSetbits()) {
                tmpReturnArray[tmpPositivePosition] = 1;
            }
        }
        return tmpReturnArray;
    }
    //
    /**
     * Returns bit array for the specified map. The map represents a molecule based on its fragments, which are
     * represented by unique SMILES in the key set and whose frequencies are mapped in the value set.
     * But the map can also contain arbitrary fragment sets. This method is a convenience method and
     * the given frequencies are not used.
     * And the  method is only available for bit fingerprints based on unique SMILES comparisons.
     * @see #getBitArray(List)
     *
     * @param aUniqueSmilesToFrequencyMap  map usually represents a molecule by representing the fragments of
     * the molecule by unique SMILES in the key set and indicating their frequency in the value set. In principle,
     * however,such a map can be applied to any set of fragments.
     * @return int[] bit array
     * @throws NullPointerException is thrown if the map aUniqueSmilesToFrequencyMap is
     * null or contains keys or values that are null respectively.
     * @throws IllegalArgumentException is thrown if the map aUniqueSmilesToFrequencyMap
     * contains keys or values that are blank/empty, respectively.
     */
    public int[] getBitArray(Map<String,Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        //converts BitSetFingerprint to BitArray
        int[] tmpReturnArray = new int[this.fragmentFingerprintMasterVectorArray.length];
        if (this.cacheFingerprintSetting && this.cacheBitSetFingerprint != null) {
            for (int tmpPositivePosition : this.cacheBitSetFingerprint.getSetbits()) {
                tmpReturnArray[tmpPositivePosition] = 1;
            }
        } else {
            List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
            for (Map.Entry<String, Integer> tmpEntry : aUniqueSmilesToFrequencyMap.entrySet()) {
                if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                    throw new NullPointerException("Given map of string and integer instances contains " +
                            "instances that are null.");
                } else if (tmpEntry.getKey().isEmpty() || tmpEntry.getKey().isBlank()) {
                    throw new IllegalArgumentException("Given map of strings an integer instances contains strings that are blank/empty.");
                }
                tmpListOfUniqueSmiles.add(tmpEntry.getKey());
            }
            for (int tmpPositivePosition : this.getBitFingerprint(tmpListOfUniqueSmiles).getSetbits()) {
                tmpReturnArray[tmpPositivePosition] = 1;
            }
        }
        return tmpReturnArray;
    }
    //
    /**
     * Returns a CountArray, which is created based on the given parameter. The map represents a molecule based on its
     * fragments, which are represented by unique SMILES in the key set and whose frequencies are mapped in
     * the value set. But the map can also contain arbitrary fragment sets.
     * The size of the array corresponds to the number of predefined (key) fragments passed during initialization.
     * However, the size may differ if there are duplicates in the specified predefined fragments, as they
     * will be ignored/removed.
     * This method is only available for count fingerprints based on unique SMILES comparisons.
     *
     * @param aUniqueSmilesToFrequencyMap map usually represents a molecule by representing the fragments of
     * the molecule by unique SMILES in the key set and indicating their frequency in the value set. In principle,
     * however,such a map can be applied to any set of fragments.
     * @return int[] count array
     * @throws NullPointerException is thrown if the map aUniqueSmilesToFrequencyMap is
     * null or contains keys or values that are null respectively.
     * @throws IllegalArgumentException is thrown if the map aUniqueSmilesToFrequencyMap
     * contains keys or values that are blank/empty, respectively.
     */
    public int[] getCountArray(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        for (Map.Entry<String, Integer> tmpEntry : aUniqueSmilesToFrequencyMap.entrySet()) {
            if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                throw new NullPointerException("Given map of string and integer instances contains " +
                        "instances that are null.");
            } else if (tmpEntry.getKey().isEmpty() || tmpEntry.getKey().isBlank()) {
                throw new IllegalArgumentException("Given map of strings an integer instances contains strings that are blank/empty.");
            }
            for(int i = 1; i <= tmpEntry.getValue(); i++) {
                tmpListOfUniqueSmiles.add(tmpEntry.getKey());
            }
        }
        return this.createCountArray(tmpListOfUniqueSmiles);
    }
    //
    /**
     * Returns the count array for the specified list.
     * This method is only available for count fingerprints based on unique SMILES comparisons.
     * @see #getCountArray(Map)
     *
     * @param aListOfUniqueSmiles is a list that stores molecule fragments or arbitrary fragments
     * in the form of unique SMILES.
     * @return int[] count array
     * @throws NullPointerException is thrown if the list aListOfUniqueSmiles is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    public int[] getCountArray(List<String> aListOfUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        return this.createCountArray(aListOfUniqueSmiles);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Generates count array for the specified list (molecule).
     * Among other things, already generated results are used to generate the array.
     * For example, if a count fingerprint has already been generated for the given list of unique SMILES or for
     * the given molecule, the result of the count fingerprint is expanded into an array. Otherwise,
     * the count fingerprint is generated first and then the count array.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * @return int[] count array
     */
    private int[] createCountArray(List<String> aListOfUniqueSmiles) {
        int[] tmpCountArray = new int[this.uniqueSmilesToPositionMap.size()];
        CountFingerprint tmpCountFingerprint;
        if (this.cacheFingerprintSetting) {
            if (this.cacheCountFingerprint == null) {
                this.cacheCountFingerprint = (CountFingerprint) this.getCountFingerprint(aListOfUniqueSmiles);
            }
            tmpCountFingerprint = this.cacheCountFingerprint;
        } else {
            tmpCountFingerprint = (CountFingerprint) this.getCountFingerprint(aListOfUniqueSmiles);
        }
        Map<Integer, Integer> tmpPositionToFrequencyMap = tmpCountFingerprint.getSmilesPositionToFrequencyMap();
        for (Map.Entry<Integer, Integer> tmpEntry : tmpPositionToFrequencyMap.entrySet()) {
            tmpCountArray[tmpEntry.getKey()] = tmpEntry.getValue();
        }
        return tmpCountArray;
    }
    //
    /**
     * Method stores all key fragments specified during initialization in an array. It ensures that there are
     * no fragment duplicates in the array.
     *
     * @return String[]
     */
    private String[] getPredefinedFragmentArrayWithoutDuplicates() {
        String[] tmpPredefinedFragmentsInArray = new String[this.uniqueSmilesToPositionMap.size()];
        for (Map.Entry<String, Integer> tmpEntry : this.uniqueSmilesToPositionMap.entrySet()) {
            tmpPredefinedFragmentsInArray[tmpEntry.getValue()] = tmpEntry.getKey();
        }
        return tmpPredefinedFragmentsInArray;
    }
    //
    /**
     * The input parameter are checked for validity.
     *
     * @param aListOfUniqueSmiles is an input list that is checked for validity.
     * @param anArgumentNullExceptionMessage NullPointerException message.
     * @param anArgumentElementNullMessage NullPointerException message for list elements.
     * @param anArgumentElementBlankEmptyMessage error message for empty/blank list elements.
     * @throws IllegalArgumentException is thrown if the input list is null.
     * @throws NullPointerException is thrown if the input list contains blank/empty strings.
     */
    private void validityCheckOfParameterList(
            List<String> aListOfUniqueSmiles,
            String anArgumentNullExceptionMessage,
            String anArgumentElementNullMessage,
            String anArgumentElementBlankEmptyMessage
    ) throws IllegalArgumentException, NullPointerException {
        Objects.requireNonNull(aListOfUniqueSmiles, anArgumentNullExceptionMessage);
        for (String tmpUniqueSmiles : aListOfUniqueSmiles) {
            Objects.requireNonNull(tmpUniqueSmiles, anArgumentElementNullMessage);
            if (tmpUniqueSmiles.isBlank()) {
                throw new IllegalArgumentException(anArgumentElementBlankEmptyMessage);
            }
        }
    }
    //
    /**
     * Private method to map a given array of fingerprint SMILES Strings to their respective positions in the fingerprint.
     * The returned HashMap includes the given SMILES as keys and the fingerprint position as values.
     *
     * @param aFragmentsArray with SMILES to build a map of
     * @return HashMap with SMILES Strings as keys and fingerprint position as value
     */
    private HashMap<String, Integer> buildUniqueSmilesToPositionMap(String[] aFragmentsArray) {
        HashMap<String, Integer> tmpUniqueSmileToPositionMap = new HashMap<>((int) (
                aFragmentsArray.length * 1.5f),
                0.75f
        );
        int tmpValuePosition = 0;
        for (String tmpKey : aFragmentsArray) {
            if (!tmpUniqueSmileToPositionMap.containsKey(tmpKey)) {
                tmpUniqueSmileToPositionMap.put(tmpKey, tmpValuePosition);
                tmpValuePosition++;
            }
        }
        return tmpUniqueSmileToPositionMap;
    }
    // </editor-fold>
}
