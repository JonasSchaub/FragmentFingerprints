/*
 * MIT License
 *
 * Copyright (c) 2026 Betuel Sevindik, Maximilian Rottmann, Felix Baensch, Jonas
 * Schaub, Christoph Steinbeck, and Achim Zielesny
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

import org.openscience.cdk.fingerprint.AbstractFingerprinter;
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
 * the form of unique SMILES (or InChIs, InChI keys, etc., alternatively) to create the fingerprint.
 * These structures must be passed as parameters when the class is instantiated (in the constructor).
 * The class implements the interface IFragmentFingerprinter,
 * which inherits the IFingerprinter (CDK).
 * The fingerprints are calculated by comparing given fragments, which are in the form of unique SMILES or any other
 * unique string representation, with the predefined fragments.
 *
 * @author Betuel Sevindik
 * @author Maximilian Rottmann
 * @author Jonas Schaub
 * @version 1.2.0.0
 */
public class FragmentFingerprinter extends AbstractFingerprinter implements IFragmentFingerprinter {
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * The fingerprint definition given during initialization is converted into a HashMap to speed up the matching
     * of the unique SMILES.
     * The Map maps the unique SMILES of the pre-defined fragments to the position they have in the fingerprint.
     */
    private final HashMap<String, Integer> uniqueSmilesToPositionMap;
    /**
     * Private integer for storing the internal map's size in order to not need on-the-fly calculation.
     */
    private final int uniqueSmilesToPositionMapSize;
    //</editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Initialization of the fragment fingerprinter by using a user-defined
     * set of fragments in the form of unique SMILES (or any unique String representation).
     * Any duplicates contained in the list will be removed, which may result in the number of predefined fragments
     * specified by the user differing from the actual number of key fragments present.
     * This means that duplicate fragment SMILES strings in the input
     * list are ignored and not represented multiple times in the fingerprint.
     *
     * @param aFragmentsForMasterVectorList in which the predefined fragments are stored.
     * @throws NullPointerException is thrown if the list param (or any of its elements) is null.
     */
    public FragmentFingerprinter(List<String> aFragmentsForMasterVectorList) throws NullPointerException {
        if (aFragmentsForMasterVectorList == null) {
            throw new NullPointerException("aFragmentsForMasterVectorList (list of string instances) is null.");
        }
        for (String fragment : aFragmentsForMasterVectorList) {
            if (fragment == null) {
                throw new NullPointerException("aFragmentsForMasterVectorList (at least one list element) is null.");
            }
        }
        this.uniqueSmilesToPositionMap = this.buildUniqueSmilesToPositionMap(aFragmentsForMasterVectorList);
        this.uniqueSmilesToPositionMapSize = this.uniqueSmilesToPositionMap.size();
    }
    /**
     * Initialization of the fragment fingerprinter by using a user-defined
     * set of fragments in the form of unique SMILES (or any unique String representation).
     * Any duplicates contained in the array will be removed, which may result in the number of predefined fragments
     * specified by the user differing from the actual number of key fragments present.
     * This means that duplicate fragment SMILES strings in the input
     * array are ignored and not represented multiple times in the fingerprint.
     *
     * @param aFragmentsForMasterVectorArray in which the predefined fragments are stored.
     * @throws NullPointerException is thrown if the array param (or any of its elements) is null.
     */
    public FragmentFingerprinter(String[] aFragmentsForMasterVectorArray) throws NullPointerException {
        if (aFragmentsForMasterVectorArray == null) {
            throw new NullPointerException("aFragmentsForMasterVectorArray (array of string instances) is null.");
        }
        for (String fragment : aFragmentsForMasterVectorArray) {
            if (fragment == null) {
                throw new NullPointerException("aFragmentsForMasterVectorArray (at least one array element) is null.");
            }
        }
        this.uniqueSmilesToPositionMap = this.buildUniqueSmilesToPositionMap(aFragmentsForMasterVectorArray);
        this.uniqueSmilesToPositionMapSize = this.uniqueSmilesToPositionMap.size();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="IFragmentFingerprinter methods">
    /**
     * Method to generate the bit fingerprint.
     * A given list of unique SMILES is compared with the predefined fragments.
     * If there is a match, the position of the unique SMILES is determined from the map and set to true in the
     * initialized BitSet. The method is intended to generate the fingerprint for one molecule but can in principle be
     * applied to any fragment set, e.g. originating from a cluster of multiple molecules.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * @return BitSetFingerprint. BitSetFingerprint is a CDK class that implements the IBitFingerprint interface of CDK.
     * @throws NullPointerException is thrown if the list param (or any of its elements) is null.
     */
    @Override
    public IBitFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles) throws NullPointerException {
        if  (aListOfUniqueSmiles == null) {
            throw new NullPointerException("Given list of string instances is null.");
        }
        for (String fragment : aListOfUniqueSmiles) {
            if (fragment == null) {
                throw new NullPointerException("Given list includes at least one null element.");
            }
        }
        BitSet tmpBitSet = new BitSet(this.uniqueSmilesToPositionMapSize);
        if (aListOfUniqueSmiles.isEmpty()) {
            return new BitSetFingerprint(tmpBitSet);
        }
        for (String tmpSmiles : aListOfUniqueSmiles) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpSmiles)) {
                tmpBitSet.set(this.uniqueSmilesToPositionMap.get(tmpSmiles), true);
            }
        }
        return new BitSetFingerprint(tmpBitSet);
    }
    //
    /**
     * Generates count fingerprint for a molecule based on its fragments represented by unique SMILES strings
     * in the key set and their frequencies in the value set of the given map. Given fragment SMILES codes that
     * are not part of the set given at initialization of this class, are ignored. The frequencies of those matching
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
     */
    @Override
    public ICountFingerprint getCountFingerprint(Map<String, Integer> aSmilesToFrequencyMap) throws NullPointerException {
        Objects.requireNonNull(aSmilesToFrequencyMap, "Given map of string and integer instances is null.");
        HashMap<Integer, Integer> tmpPositionToFrequencyMap = new HashMap<>(
                (int) (this.uniqueSmilesToPositionMapSize / 0.75 + 1),
                0.75f);
        for (Map.Entry<String, Integer> tmpEntry : aSmilesToFrequencyMap.entrySet()) {
            if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                throw new NullPointerException("Given map of string and integer instances contains " +
                        "instances that are null.");
            } else if (this.uniqueSmilesToPositionMap.containsKey(tmpEntry.getKey())) {
                //easier debugging (speaking from experience)
                int tmpPosition = this.uniqueSmilesToPositionMap.get(tmpEntry.getKey());
                int tmpFrequency = tmpEntry.getValue();
                tmpPositionToFrequencyMap.put(tmpPosition, tmpFrequency);
            }
        }
        return new CountFingerprint(this.uniqueSmilesToPositionMapSize, tmpPositionToFrequencyMap);
    }
    //
    /**
     * Generates a count fingerprint for a molecule based on its fragments, represented by unique SMILES in the
     * list given as parameters. Given fragment SMILES codes that are not part of the set given at initialization of
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
     * @throws NullPointerException is thrown if the list param (or any of its elements) is null.
     */
    @Override
    public ICountFingerprint getCountFingerprint(List<String> aUniqueSmilesList) throws NullPointerException {
        HashMap<String, Integer> tmpUniqueSmilesToFrequencyCountMap
                = new HashMap<>((int) (this.uniqueSmilesToPositionMapSize / 0.75 + 1), 0.75f);
        Objects.requireNonNull(aUniqueSmilesList, "aUniqueSmilesToFrequencyList (list of string instances) is null.");
        for (String tmpSmiles : aUniqueSmilesList) {
            Objects.requireNonNull(tmpSmiles, "aUniqueSmilesToFrequencyList (at least one list element) is null.");
            if (!tmpUniqueSmilesToFrequencyCountMap.containsKey(tmpSmiles)) {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, 1);
            } else {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, tmpUniqueSmilesToFrequencyCountMap.get(tmpSmiles) + 1);
            }
        }
        return this.getCountFingerprint(tmpUniqueSmilesToFrequencyCountMap);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="IFingerprinter methods">
    /**
     * {@inheritDoc}
     *
     * @throws UnsupportedOperationException as the method is no longer supported.
     */
    @Override
    public BitSet getFingerprint(IAtomContainer mol) throws UnsupportedOperationException {
        throw new UnsupportedOperationException("This method is no longer supported. " +
                "Please use the CDK SubstructureFingerprinter class for substructure-based fingerprinting instead.");
    }
    //
    /**
     * {@inheritDoc}
     * @see SubstructureFingerprinter
     *
     * @throws UnsupportedOperationException as the method is no longer supported.
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws UnsupportedOperationException {
        throw new UnsupportedOperationException("This method is no longer supported. " +
                "Please use the CDK SubstructureFingerprinter class for substructure-based fingerprinting instead.");
    }
    //
    /**
     * {@inheritDoc}
     * @see SubstructureFingerprinter
     *
     * @throws UnsupportedOperationException as the method is no longer supported.
     */
    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws UnsupportedOperationException {
        throw new UnsupportedOperationException("This method is no longer supported. " +
                "Please use the CDK SubstructureFingerprinter class for substructure-based fingerprinting instead.");
    }
    //
    /**
     * UnsupportedOperationException. This method is not supported.
     * {@inheritDoc}
     *
     * @throws UnsupportedOperationException method is not supported
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws UnsupportedOperationException {
       throw new UnsupportedOperationException("This method is no longer supported.");
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since the FragmentFingerprinter is a key-based fingerprint, the size of the fingerprint is equal
     * to the number of predefined fragments (unique SMILES) passed at initialization (minus removed duplicates).
     *
     * @return int
     */
    @Override
    public int getSize() {
        return this.uniqueSmilesToPositionMapSize;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Returns the bit definitions i.e. which bit stands for which fragment SMILES.
     * Important, the number of possible bit definitions may differ from the number of key
     * fragments passed during initialization, since duplicates are removed.
     * Equally important, be advised that this operation runs in O(n) time.
     *
     * @param aBitPosition in the fingerprint.
     * @return unique SMILES corresponding to the specified position.
     * @throws IllegalArgumentException is thrown if the given bit position is not (or cannot be) present in the fingerprint.
     */
    public String getBitDefinition(int aBitPosition) throws IllegalArgumentException {
        if (aBitPosition < 0) {
            throw new IllegalArgumentException("Given bit position cannot be smaller than 0.");
        } else if (aBitPosition >= this.uniqueSmilesToPositionMapSize) {
            throw new IllegalArgumentException("Given bit position cannot be larger than position map.");
        } else {
            for (Map.Entry<String, Integer> tmpEntry : this.uniqueSmilesToPositionMap.entrySet()) {
                if (tmpEntry.getValue() == aBitPosition) {
                    return tmpEntry.getKey();
                }
            }
        }
        throw new IllegalArgumentException("This bit is not defined/present in the fingerprint.");
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
     * @throws NullPointerException is thrown if the list param (or any of its elements) is null.
     */
    public int[] getBitArray(List<String> aListOfUniqueSmiles) throws NullPointerException {
        if (aListOfUniqueSmiles == null) {
            throw new NullPointerException("aListOfUniqueSmiles (list of string instances) is null.");
        }
        for (String tmpSmiles : aListOfUniqueSmiles) {
            if (tmpSmiles == null) {
                throw new NullPointerException("aListOfUniqueSmiles (at least one list element) is null.");
            }
        }
        //converts BitSetFingerprint to BitSet
        int[] tmpReturnArray = new int[this.uniqueSmilesToPositionMapSize];
        for (int tmpPositivePosition : this.getBitFingerprint(aListOfUniqueSmiles).getSetbits()) {
            tmpReturnArray[tmpPositivePosition] = 1;
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
     */
    public int[] getBitArray(Map<String,Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        //converts BitSetFingerprint to BitArray
        int[] tmpReturnArray = new int[this.uniqueSmilesToPositionMapSize];
        for (String tmpKey : aUniqueSmilesToFrequencyMap.keySet()) {
            if (tmpKey == null) {
                throw new NullPointerException("Map contains null key.");
            }
            Integer tmpPosition = this.uniqueSmilesToPositionMap.get(tmpKey);
            if (tmpPosition != null) {
                tmpReturnArray[tmpPosition] = 1;
            }
        }
        return tmpReturnArray;
    }
    //
    /**
     * Method directly returning the BitSet of the fingerprint generated from the given list of SMILES.
     *
     * @param aListOfUniqueSmiles storing the fragments or molecules from which a fingerprint (and BitSet) is to be generated
     * @return BitSet of the generated fingerprint
     */
    public BitSet getBitSet(List<String> aListOfUniqueSmiles) throws NullPointerException {
        if (aListOfUniqueSmiles == null) {
            throw new NullPointerException("aListOfUniqueSmiles (list of string instances) is null.");
        }
        for (String tmpSmiles : aListOfUniqueSmiles) {
            if (tmpSmiles == null) {
                throw new NullPointerException("aListOfUniqueSmiles (at least one list element) is null.");
            }
        }
        return this.getBitFingerprint(aListOfUniqueSmiles).asBitSet();
    }
    /**
     * Method directly returning the BitSet of the fingerprint generated from the given frequency map.
     *
     * @param aUniqueSmilesToFrequencyMap with SMILES to frequency representation from which a fingerprint (and BitSet) is to be generated
     * @return BitSet of the generated fingerprint
     */
    public BitSet getBitSet(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap,
                "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        for (Map.Entry<String, Integer> tmpEntry : aUniqueSmilesToFrequencyMap.entrySet()) {
            if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                throw new NullPointerException("Given map of string and integer instances contains " +
                        "instances that are null.");
            }
            tmpListOfUniqueSmiles.add(tmpEntry.getKey());
        }
        return this.getBitFingerprint(tmpListOfUniqueSmiles).asBitSet();
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
     */
    public int[] getCountArray(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        int[] tmpCountArray = new int[this.uniqueSmilesToPositionMapSize];
        for (Map.Entry<String, Integer> tmpEntry : aUniqueSmilesToFrequencyMap.entrySet()) {
            if (tmpEntry.getKey() == null || tmpEntry.getValue() == null) {
                throw new NullPointerException("Given map of string and integer instances contains " +
                        "instances that are null.");
            }
            Integer tmpPosition = this.uniqueSmilesToPositionMap.get(tmpEntry.getKey());
            if (tmpPosition != null) {
                tmpCountArray[tmpPosition] = tmpEntry.getValue();
            }
        }
        return tmpCountArray;
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
     * @throws NullPointerException is thrown if the list param (or any of its elements) is null.
     */
    public int[] getCountArray(List<String> aListOfUniqueSmiles) throws NullPointerException {
        if (aListOfUniqueSmiles == null) {
            throw new NullPointerException("aListOfUniqueSmiles (list of string instances) is null.");
        }
        for (String tmpSmiles : aListOfUniqueSmiles) {
            if (tmpSmiles == null) {
                throw new NullPointerException("aListOfUniqueSmiles (at least one list element) is null.");
            }
        }
        return this.createCountArray(aListOfUniqueSmiles);
    }
    //
    /**
     * Method to return the count/occurrences/frequency of a given SMILES String in a given CountFingerprint instance.
     *
     * <p><b>Important:</b> The {@code CountFingerprint} instance <em>must</em> have been generated by this
     * very {@code FragmentFingerprinter} instance (or by a {@code FragmentFingerprinter} with an identical
     * configuration and fragment/SMILES index mapping). The method assumes that the internal mapping from
     * unique SMILES to array positions in this fingerprinter matches the mapping used when the given
     * {@code CountFingerprint} was created.</p>
     * <p>If this requirement is violated, the lookup may return incorrect counts for the requested SMILES or,
     * if the position is not present in the fingerprint, may result in an {@link IllegalArgumentException}
     * or other runtime errors. There is no automatic runtime compatibility check; callers are responsible
     * for ensuring that only compatible {@code CountFingerprint} instances are passed in.</p>
     * <p><b>Correct usage:</b> Always obtain the {@code CountFingerprint} to be queried from methods of the
     * same {@code FragmentFingerprinter} instance that is used to call this method (or from another instance
     * that is guaranteed to use the exact same fragment set and SMILES-to-position mapping).</p>
     *
     * @param aSmiles String to get count for
     * @param aCountFingerprint wherein to search for SMILES String
     * @return integer count of given SMILES String
     * @throws IllegalArgumentException if SMILES is not present in fingerprint
     * @throws NullPointerException if either/both params are null
     */
    public int count(String aSmiles, CountFingerprint aCountFingerprint) throws IllegalArgumentException {
        Objects.requireNonNull(aSmiles, "Given SMILES was null.");
        Objects.requireNonNull(aCountFingerprint, "Given CountFingerprint was null.");
        if(!this.uniqueSmilesToPositionMap.containsKey(aSmiles)) {
            throw new IllegalArgumentException("The given SMILES string is not present in defined fingerprint.");
        }
        int tmpPosition =  this.uniqueSmilesToPositionMap.get(aSmiles);
        try {
            return aCountFingerprint.getSmilesPositionToFrequencyMap().get(tmpPosition);
        } catch (NullPointerException aNullpointerException) {
            throw new IllegalArgumentException("Given SMILES string does not occur in given CountFingerprint.");
        }
    }
    //
    /**
     * Public method to get a float[][] matrix containing the fingerprints for the given maps list (fragment sets
     *  of distinct molecules) in regard to the pre-defined fragment fingerprint.
     * It shall be noted that each Map in the given List represents ONE molecule's fragment frequencies (fragment-frequency pairs).
     * Therefore, the whole list represents the entirety of fragments.
     * Further, a setting whether to use "bit set" or "count/frequency of fragment" is available.
     * It is recommended to check the size of the pre-initialized float matrix. In case of a matrix with a pre-initialized
     *  size smaller than the actual required size, an exception will be thrown and no matrix filling will take place.
     *
     * @param aFragmentsFrequenciesMapsArray contains maps of SMILES-Frequency pairs of fragments
     * @param aFloatDataMatrix is a pre-initialized float[][] matrix to write data into
     * @param anUseBitArrayStatement is the setting whether bit or count/frequency representation should be used
     *
     * @throws NullPointerException if given array of frequencies maps is null
     * @throws IllegalArgumentException if matrix sizes (row or column) are insufficient
     */
    public void getFragmentsComponentsFloatMatrix(
            //array of maps of all molecules' fragments of a specific fragmentation
            Map<String, Integer>[] aFragmentsFrequenciesMapsArray,
            float[][] aFloatDataMatrix,
            boolean anUseBitArrayStatement) {
        Objects.requireNonNull(aFragmentsFrequenciesMapsArray);
        //checks whole matrix (row count) for allowed size by comparing to given list size
        if (aFragmentsFrequenciesMapsArray.length > aFloatDataMatrix.length) {
            throw new IllegalArgumentException("Given fragments maps array length was larger than provided matrix' length.");
        } else {
            //checks each matrix array (~column count) for allowed size
            for (float[] tmpInnerFloatArray : aFloatDataMatrix) {
                if (this.uniqueSmilesToPositionMapSize > tmpInnerFloatArray.length) {
                    throw new IllegalArgumentException("Given matrix' array (row) size was smaller than the pre-defined fingerprint.");
                }
            }
            //float[row count][column count]
            for (int i = 0; i < aFragmentsFrequenciesMapsArray.length; i++) {
                this.getFloatFingerprint(
                        aFragmentsFrequenciesMapsArray[i], //gives frequencies map of specific molecule
                        aFloatDataMatrix[i], //defines where to put generated fingerprint in the matrix
                        anUseBitArrayStatement); //defines whether bit or count/frequency representation is used
            }
        }
    }
    //</editor-fold>
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
        int[] tmpCountArray = new int[this.uniqueSmilesToPositionMapSize];
        CountFingerprint tmpCountFingerprint = (CountFingerprint) this.getCountFingerprint(aListOfUniqueSmiles);
        Map<Integer, Integer> tmpPositionToFrequencyMap = tmpCountFingerprint.getSmilesPositionToFrequencyMap();
        for (Map.Entry<Integer, Integer> tmpEntry : tmpPositionToFrequencyMap.entrySet()) {
            tmpCountArray[tmpEntry.getKey()] = tmpEntry.getValue();
        }
        return tmpCountArray;
    }
    //
    /**
     * Private method to map a given array of fingerprint SMILES Strings to their respective positions in the fingerprint.
     * The returned HashMap includes the given SMILES as keys and the fingerprint position as values.
     *
     * @param aFragmentsList with SMILES to build a map of
     * @return HashMap with SMILES Strings as keys and fingerprint position as value
     */
    private HashMap<String, Integer> buildUniqueSmilesToPositionMap(List<String> aFragmentsList) {
        HashMap<String, Integer> tmpUniqueSmileToPositionMap = new HashMap<>(
                (int) (aFragmentsList.size() / 0.75 + 1), 0.75f);
        int tmpValuePosition = 0;
        for (String tmpString : aFragmentsList) {
            if (!tmpUniqueSmileToPositionMap.containsKey(tmpString)) {
                tmpUniqueSmileToPositionMap.put(tmpString, tmpValuePosition);
                tmpValuePosition++;
            }
        }
        return tmpUniqueSmileToPositionMap;
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
                aFragmentsArray.length / 0.75 + 1),0.75f);
        int tmpValuePosition = 0;
        for (String tmpString : aFragmentsArray) {
            if (!tmpUniqueSmileToPositionMap.containsKey(tmpString)) {
                tmpUniqueSmileToPositionMap.put(tmpString, tmpValuePosition);
                tmpValuePosition++;
            }
        }
        return tmpUniqueSmileToPositionMap;
    }
    //
    /**
     * Method to generate the float fingerprint of a given molecule's fragments (i.e. the given map).
     * A setting is defined whether bit or count/frequency representation shall be used for generating the fingerprint.
     *
     * @param aFragmentsFrequenciesMap contains SMILES-frequency pairs of the molecule's fragments
     * @param aPreInitFloatArray is a pre-initialized float[] array (-> matrix row)
     * @param anUseBitArrayStatement is the setting whether bit or count/frequency representation is to be used
     *
     * @throws NullPointerException if parameter fragments map is null
     * @throws IllegalArgumentException if the frequency of a fragment was smaller than 0
     */
    private void getFloatFingerprint(
            Map<String, Integer> aFragmentsFrequenciesMap,
            float[] aPreInitFloatArray,
            boolean anUseBitArrayStatement
    ) {
        Objects.requireNonNull(aFragmentsFrequenciesMap, "Given map was null.");
        //validity check of frequencies
        for (int tmpFrequency : aFragmentsFrequenciesMap.values()) {
            if (tmpFrequency < 0) {
                throw new IllegalArgumentException("Frequency of a given fragment was smaller than 0.");
            }
        }
        for (Map.Entry<String, Integer> tmpFragmentEntry : aFragmentsFrequenciesMap.entrySet()) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpFragmentEntry.getKey())) {
                if (anUseBitArrayStatement) {
                    aPreInitFloatArray[this.uniqueSmilesToPositionMap.get(tmpFragmentEntry.getKey())] = 1.0f;
                } else {
                    aPreInitFloatArray[this.uniqueSmilesToPositionMap.get(tmpFragmentEntry.getKey())] = tmpFragmentEntry.getValue();
                }
            }
        }
    }
    // </editor-fold>
}
