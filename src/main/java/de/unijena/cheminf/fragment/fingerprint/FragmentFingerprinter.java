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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Class to generate fragment fingerprints. Bit and count fragment fingerprints can be generated.
 * Fragment fingerprints are key-based fingerprints.
 * Thus, the class requires predefined structures/fragments in
 * the form of unique SMILES to create the fingerprint. These structures must be passed when the
 * class is instantiated (in the constructor). The class implements the interface IFragmentFingerprinter,
 * which inherits the IFingerprinter (CDK), which allows the class to compute fingerprints in 2 ways.
 * The first way to calculate a bit or count fingerprint is to perform a substructure comparison with all
 * predefined fragments for a given IAtomContainer. The fingerprint created by the substructure search is based on
 * the CDK class SubtructureFingerprinter. The predefined fragment SMILES are interpreted as SMARTS patterns by the
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
    private final String[] fragmentsForBitSetArray;
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
     * The fragmentArray is converted into a HashMap to speed up the matching of the unique SMILES.
     * The Map maps the unique SMILES of the predefined fragments to the position they have in the array.
     */
    private HashMap<String, Integer> uniqueSmilesToPositionMap;
    /**
     * Bit fingerprint for storing the calculated fragment bit fingerprint.
     */
    private BitSetFingerprint cacheBitSetFingerprint;
    /**
     * This map is a raw map. If a match between key fragments and molecule fragments (or any fragments) occurs
     * during the creation of the count fingerprint,the position of the key fragment in the fingerprint is mapped to
     * the frequency of occurrence of the fragment.
     */
    private HashMap<Integer,Integer> cacheRawCountMap;
    /**
     * The list is a clone of the list passed as a parameter when creating the bit fingerprint.
     * It is used to check whether a bit fingerprint already exists for a given list of fragments or molecule fragments.
     */
    private List<String> cacheListToGenerateBitFingerprint;
    /**
     * The list stores all keys of a map that are necessary to create the count fingerprint.
     * The map usually represents a molecule by representing the fragments of the molecule by unique SMILES in the
     * key set and indicating their frequency in the value set. In principle, however,such a map can be applied to any
     * set of fragments.This list is used to check whether a count fingerprint has already been created for the
     * specified map.
     */
    private List<String> cacheListToGenerateCountFingerprint;
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
     * @param aFragmentForBitSetList is the ist in which the predefined fragments are stored.
     * @throws NullPointerException is thrown if the list aFragmentForBitSetList is null.
     * @throws IllegalArgumentException is thrown if the list contains blank strings.
     */
    //ToDo: class variables are not good here
    public FragmentFingerprinter(List<String> aFragmentForBitSetList) throws NullPointerException, IllegalArgumentException {
        // Check whether aFragmentForBitSetList is null or whether there are elements (strings) in the list that are empty.
        this.validityCheckOfParameterList(aFragmentForBitSetList,"aFragmentForBitSetList (list of string instances) is null.",
                "aFragmentForBitSetList (at least one list element) is null.",
                "aFragmentForBitSetList (at least one list element) is blank/empty.");
        this.fragmentsForBitSetArray = aFragmentForBitSetList.toArray(new String[aFragmentForBitSetList.size()]);
        this.uniqueSmilesToPositionMap = this.buildUniqueSmilesToPositionMap(this.fragmentsForBitSetArray);
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
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        BitSet tmpBitSet = new BitSet(this.uniqueSmilesToPositionMap.size());
        Set<String> tmpUniqueSmilesSet = new HashSet<>((int) (aListOfUniqueSmiles.size() * 1.5f));
        tmpUniqueSmilesSet.addAll(aListOfUniqueSmiles);
        ArrayList<String> tmpUniqueSmilesList = new ArrayList<>(tmpUniqueSmilesSet);
        this.cacheListToGenerateBitFingerprint = (ArrayList<String>) tmpUniqueSmilesList.clone();
        for(String tmpUniqueSmilesWithoutDuplicates : tmpUniqueSmilesList) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpUniqueSmilesWithoutDuplicates)) {
                int tmpPosition = uniqueSmilesToPositionMap.get(tmpUniqueSmilesWithoutDuplicates);
                tmpBitSet.set(tmpPosition, true);
            }
        }
        this.cacheBitSetFingerprint = new BitSetFingerprint(tmpBitSet);
        return this.cacheBitSetFingerprint;
    }

    //ToDo: implicitly with method below -> setFingerprint() necessary
    public BitSetFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles, Map<String, Integer> aSmilesToPositionMap) {
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        BitSet tmpBitSet = new BitSet(aSmilesToPositionMap.size());
        Set<String> tmpUniqueSmilesSet = new HashSet<>((int) (aListOfUniqueSmiles.size() * 1.5f));
        //hash set is used to deduplicate
        tmpUniqueSmilesSet.addAll(aListOfUniqueSmiles);
        for (String tmpSmiles : tmpUniqueSmilesSet) {
            if (aSmilesToPositionMap.containsKey(tmpSmiles)) {
                int tmpPosition = aSmilesToPositionMap.get(tmpSmiles);
                tmpBitSet.set(tmpPosition, true);
            }
        }
        return new BitSetFingerprint(tmpBitSet);
    }

    public float[] getFloatBitFingerprint(
            //aFragmentsUniqueSmilesList contains fragments of ONE molecule
            List<String> aFragmentsUniqueSmilesList,
            Map<String, Integer> aSmilesToPositionMap,
            float[] aPreInitFloatArray)
    {
        //ToDo: validity check necessary at this point? especially when partially "empty" arrays are passed through?
//        this.validityCheckOfParameterList(aFragmentsUniqueSmilesList,"aFragmentsUniqueSmilesList (list of string instances) is null.",
//                "aFragmentsUniqueSmilesList (at least one list element) is null.",
//                "aFragmentsUniqueSmilesList (at least one list element) is blank/empty.");
        Set<String> tmpUniqueSmilesSet = new HashSet<>((int) (aFragmentsUniqueSmilesList.size() * 1.5f));
        tmpUniqueSmilesSet.addAll(aFragmentsUniqueSmilesList);
        for (String tmpSmiles : tmpUniqueSmilesSet) {
            if (aSmilesToPositionMap.containsKey(tmpSmiles)) {
                aPreInitFloatArray[aSmilesToPositionMap.get(tmpSmiles)] = 1.0f;
            }
        }
        return aPreInitFloatArray;
    }

    public float[] getFloatCountFingerprint(
            List<String> aFragmentsUniqueSmilesList,
            Map<String, Integer> aSmilesToPositionMap,
            float[] aPreInitFloatArray)
    {
        this.validityCheckOfParameterList(aFragmentsUniqueSmilesList,"aListOfFragmentsUniqueSmiles (list of string instances) is null.",
                "aListOfFragmentsUniqueSmiles (at least one list element) is null.",
                "aListOfFragmentsUniqueSmiles (at least one list element) is blank/empty.");
        for (String tmpSmiles : aFragmentsUniqueSmilesList) {
            if (aSmilesToPositionMap.containsKey(tmpSmiles)) {
                aPreInitFloatArray[aSmilesToPositionMap.get(tmpSmiles)]++;
            }
        }
        return aPreInitFloatArray;
    }
    //
    /**
     * Generates count fingerprint for a molecule based on its fragments represented by unique SMILES strings
     * in the key set and their frequencies in the value set of the given map. Given fragment SMILES codes that
     * are not part of the set given at initialisation of this class, are ignored. The frequencies of those matching
     * with the predefined set are used to construct the fingerprint. The method is intended to generate the fingerprint
     * for one molecule but can in principle be applied to any fragment set, e.g. originating from a cluster of
     * multiple molecules
     *
     * @param aUniqueSmilesToFrequencyMap map usually represents a molecule by representing the fragments of
     * the molecule by unique SMILES in the key set and indicating their frequency in the value set. In principle,
     * however,such a map can be applied to any set of fragments.
     * To be able to calculate the fingerprint for a molecule, the fragments must belong to a molecule.
     * @return count fingerprint
     * @throws NullPointerException  is thrown if the map aUniqueSmilesToFrequencyMap is
     * null or contains keys or values that are null respectively.
     * @throws IllegalArgumentException is thrown if the map aUniqueSmilesToFrequencyMap
     * contains keys or values that are blank/empty, respectively.
     */
    @Override
    public CountFingerprint getCountFingerprint(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException,IllegalArgumentException {
        this.cacheRawCountMap = new HashMap<>((int) (this.uniqueSmilesToPositionMap.size() * 1.5f), 0.75f);
        this.cacheListToGenerateCountFingerprint = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        for (String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpUniqueSmiles == null || aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null.");
            }
            if(tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyMap (Map of strings an integer instances) contains strings that are blank/empty.");
            }
            if (this.uniqueSmilesToPositionMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = this.uniqueSmilesToPositionMap.get(tmpUniqueSmiles);
                this.cacheRawCountMap.put(tmpPosition,aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles));
            }
            for(int i = 1; i<=aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles); i++) {
                this.cacheListToGenerateCountFingerprint.add(tmpUniqueSmiles);
            }
        }
        return new CountFingerprint(this.fragmentsForBitSetArray, this.cacheRawCountMap);
    }
    //
    /**
     * Generates a count fingerprint for a molecule based on its fragments, represented by unique SMILES in the
     * list given as parameters. Given fragment SMILES codes that are not part of the set given at initialisation of
     * this class, are ignored. The frequencies of those matching with the predefined set are used to construct the
     * fingerprint The frequency of individual fragments depends on how often they occur in the specified list.
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
            if(tmpSmiles.isBlank() || tmpSmiles.isEmpty()) {
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
        StringBuilder tmpVersionDescriptionBuilder = new StringBuilder();
        tmpVersionDescriptionBuilder.append(getClass().getSimpleName()).append("/").append(FragmentFingerprinter.FRAGMENT_FINGERPRINTER_VERSION)
                .append(' ').append("num_bits").append('=').append(this.uniqueSmilesToPositionMap.size());
        return tmpVersionDescriptionBuilder.toString();
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
        Set<String> tmpUniqueSmilesSet = new HashSet<>((int) (aListOfUniqueSmiles.size() * 1.5f));
        tmpUniqueSmilesSet.addAll(aListOfUniqueSmiles);
        List<String> tmpListWithoutPossibleDuplicates = new ArrayList<>(tmpUniqueSmilesSet);
        return this.createBitArray(tmpListWithoutPossibleDuplicates);
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
        List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        for(String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpUniqueSmiles == null || aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null.");
            }
            if(tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyMap (Map of strings an integer instances) contains strings that are blank/empty.");
            }
            tmpListOfUniqueSmiles.add(tmpUniqueSmiles);
        }
        return this.createBitArray(tmpListOfUniqueSmiles);
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
        for(String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpUniqueSmiles == null || aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null.");
            }
            if(tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyMap (Map of strings an integer instances) contains strings that are blank/empty.");
            }
            for(int i = 1; i<=aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles); i++) {
                tmpListOfUniqueSmiles.add(tmpUniqueSmiles);
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
     * @return float[] bit array
     * @throws NullPointerException is thrown if the map aUniqueSmilesToFrequencyMap is
     * null or contains keys or values that are null respectively.
     * @throws IllegalArgumentException is thrown if the map aUniqueSmilesToFrequencyMap
     * contains keys or values that are blank/empty, respectively.
     */
    public float[] getBitFloatArray(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        for(String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpUniqueSmiles == null || aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null.");
            }
            if(tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyMap (Map of strings an integer instances) contains strings that are blank/empty.");
            }
            tmpListOfUniqueSmiles.add(tmpUniqueSmiles);
        }
        return this.createFloatBitArray(tmpListOfUniqueSmiles, this.cacheBitSetFingerprint);
    }
    /**
     * Returns bit array for specified list.
     * The size of the array corresponds to the number of predefined (key) fragments passed during initialization.
     * However, the size may differ if there are duplicates in the specified predefined fragments, as they
     * will be ignored/removed.
     * This method is only available for bit fingerprints based on unique SMILES comparisons.
     *
     * @param aListOfUniqueSmiles is a list that stores molecule fragments or arbitrary fragments
     * in the form of unique SMILES.
     * @return float[] bit array
     * @throws NullPointerException is thrown if the list aListOfUniqueSmiles is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    public float[] getBitFloatArray(List<String> aListOfUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        this.validityCheckOfParameterList(aListOfUniqueSmiles,"aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        Set<String> tmpUniqueSmilesSet = new HashSet<>((int) (aListOfUniqueSmiles.size() * 1.5f));
        tmpUniqueSmilesSet.addAll(aListOfUniqueSmiles);
        List<String> tmpListWithoutPossibleDuplicates = new ArrayList<>(tmpUniqueSmilesSet);
        return this.createFloatBitArray(tmpListWithoutPossibleDuplicates, this.cacheBitSetFingerprint);
    }
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
    public float[] getCountFloatArray(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        List<String> tmpListOfUniqueSmiles = new ArrayList<>(aUniqueSmilesToFrequencyMap.size());
        for(String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpUniqueSmiles == null || aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null.");
            }
            if(tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException("aUniqueSmilesToFrequencyMap (Map of strings an integer instances) contains strings that are blank/empty.");
            }
            for(int i = 1; i<=aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles); i++) {
                tmpListOfUniqueSmiles.add(tmpUniqueSmiles);
            }
        }
        return this.createCountFloatArray(tmpListOfUniqueSmiles);
    }
    /**
     * Returns the count array for the specified list.
     * This method is only available for count fingerprints based on unique SMILES comparisons.
     * @see #getCountArray(Map)
     *
     * @param aListOfUniqueSmiles is a list that stores molecule fragments or arbitrary fragments
     * in the form of unique SMILES.
     * @return float[] count array
     * @throws NullPointerException is thrown if the list aListOfUniqueSmiles is null.
     * @throws IllegalArgumentException is thrown if the list aListOfUniqueSmiles contains blank/empty strings.
     */
    public float[] getCountFloatArray(List<String> aListOfUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        this.validityCheckOfParameterList(aListOfUniqueSmiles, "aListOfUniqueSmiles (list of string instances) is null.",
                "aListOfUniqueSmiles (at least one list element) is null.",
                "aListOfUniqueSmiles (at least one list element) is blank/empty.");
        return this.createCountFloatArray(aListOfUniqueSmiles);
    }
    /**
     * Returns the float matrix for the specified list of molecules.
     * A list of fragments must be given to generate the fingerprints.
     * Additional components for the fingerprint may be specified (added to the generated fingerprint array). These must
     * be in the same order as the given list of starter molecules!
     *
     * @param aNumberOfAdditionalComponents is an array of floats, containing the values of the generated fingerprints and descriptor values.
     * @param anUseBitArrayStatement whether the bit or count array should be used
     * @return float[][] matrix
     */
    @Deprecated
    public float[][] generateFragmentsComponentsFloatMatrix(
            //only count necessary
            HashMap<String, List<String>> anOriginToFragmentsMap,
            //count of descriptor additional columns
            float[] aNumberOfAdditionalComponents,
            boolean anUseBitArrayStatement
    ) {
        if (anUseBitArrayStatement) {
            //for matrix needed:
            // Smiles of mols = anOriginMoleculeUniqueSmilesList
            // bit array for each mol
            //for bit array needed:
            // unique smiles of fragments (without duplicates)

            //hash set is used to deduplicate, contains ALL (unique) fragments
            //the "* 5" here is completely arbitrary as no estimation for the size of fragments list can be made pre-iteration
            //ToDo: should responsibility for size lie with user? -> parameter of method?
            Set<String> tmpDeduplFragmentsSmilesSet = new HashSet<>(anOriginToFragmentsMap.size() * 5, 0.75f);
            for (Map.Entry<String, List<String>> tmpEntry : anOriginToFragmentsMap.entrySet()) {
                tmpDeduplFragmentsSmilesSet.addAll(tmpEntry.getValue());
            }
            Map<String, Integer> tmpFragmentsPositionMap = new HashMap<>((int) (tmpDeduplFragmentsSmilesSet.size() * 1.5f), 0.75f);
            int i = 0;
            for (String tmpFragmentSmiles : tmpDeduplFragmentsSmilesSet) {
                tmpFragmentsPositionMap.put(tmpFragmentSmiles, i++);
            }
            //Bit array (not frequency!)
            BitSetFingerprint[] tmpBitSetArray = new BitSetFingerprint[anOriginToFragmentsMap.size()];
            //see to do below for info on floats
            //[columns][rows]
            float[][] tmpFloatMatrix = new float[anOriginToFragmentsMap.size()][tmpDeduplFragmentsSmilesSet.size()];
            int tmpCurrentRowIndex = 0;
            //
            int j = 0;
            for (Map.Entry<String, List<String>> tmpEntry : anOriginToFragmentsMap.entrySet()) {
                String tmpOriginSmiles = tmpEntry.getKey();
                //set with fragments smiles of only ONE original molecule
                Set<String> tmpDeduplSet = new HashSet<>((int) (tmpEntry.getValue().size() * 1.5f), 0.75f);
                tmpDeduplSet.addAll(tmpEntry.getValue());
                //ToDo: Bit set and -Fingerprint exchangeable with float implementation?
                BitSet tmpBitSet = new BitSet(tmpDeduplFragmentsSmilesSet.size());
                //
                float[] tmpFragmentsFloatBitArray = new float[tmpDeduplFragmentsSmilesSet.size()];
                //
                for (String tmpFragmentSmiles : tmpDeduplSet) {
                    if (tmpFragmentsPositionMap.containsKey(tmpFragmentSmiles)) {
                        tmpBitSet.set(tmpFragmentsPositionMap.get(tmpFragmentSmiles), true);
                        //
                        tmpFragmentsFloatBitArray[tmpFragmentsPositionMap.get(tmpFragmentSmiles)] = 1.0f;
                        //
                    }
                }
                BitSetFingerprint tmpBitSetFingerprint = new BitSetFingerprint(tmpBitSet);
                tmpBitSetArray[j++] = tmpBitSetFingerprint;
                //
                //set tmpFragmentsFloatBitArray as row of tmpFloatMatrix
                //seems no other way than to iterate through whole matrix and add row
                //row index can be kept track off though, so bit quicker: no whole iteration needed -> "jump" to needed row
            }
            //result: array with BitSetFingerprints where position in array corresponds to starter molecule
            // but it's bit sets, not float arrays...
        }
        //return this.createFragmentsComponentsFloatMatrix(tmpList, aNumberOfAdditionalComponents);
        return null;
    }
    @Deprecated
    public float[][] getFragmentsComponentsFloatMatrix(
            ArrayList<List<String>> aFragmentsUniqueSmilesListsArrayList,
            int aColumnLength,
            int aRowLength,
            boolean anUseBitArrayStatement
    ) {

        return this.createFragmentsComponentsFloatMatrix(aFragmentsUniqueSmilesListsArrayList, aColumnLength, aRowLength, anUseBitArrayStatement);
    }

    public float[][] getFragmentsComponentsFloatMatrix(
            List<List<String>> aFragmentsUniqueSmilesListsArrayList,
            float[][] aFloatDataMatrix,
            boolean anUseBitArrayStatement
    ) {
        return this.fillFragmentsComponentsFloatMatrix(
                (ArrayList<List<String>>) aFragmentsUniqueSmilesListsArrayList,
                aFloatDataMatrix,
                anUseBitArrayStatement
        );
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
        if(this.cacheRawCountMap != null && this.cacheListToGenerateCountFingerprint.size() == aListOfUniqueSmiles.size()) {
            Collections.sort(aListOfUniqueSmiles);
            Collections.sort(this.cacheListToGenerateCountFingerprint);
            if(aListOfUniqueSmiles.equals(this.cacheListToGenerateCountFingerprint)) {
                for (int tmpPositivePositions : this.cacheRawCountMap.keySet()) {
                    tmpCountArray[tmpPositivePositions] = this.cacheRawCountMap.get(tmpPositivePositions);
                }
            }
        } else {
            this.cacheRawCountMap = null;
            this.cacheListToGenerateCountFingerprint = null;
            this.getCountFingerprint(aListOfUniqueSmiles);
            for (int tmpPositivePositions : this.cacheRawCountMap.keySet()) {
                tmpCountArray[tmpPositivePositions] = this.cacheRawCountMap.get(tmpPositivePositions);
            }
        }
        return tmpCountArray;
    }
    //
    /**
     * Generates bit array for the specified list (molecule).
     * Among other things, already generated results are used to generate the array.
     * For example, if a bit fingerprint has already been generated for the given list of unique SMILES or for
     * the given molecule, the result of the bit fingerprint is expanded into an array. Otherwise,
     * the bit fingerprint is generated first and then the bit array.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * @return int[] bit array
     */
    private int[] createBitArray(List<String> aListOfUniqueSmiles) {
        int[] tmpBitArray = new int[this.uniqueSmilesToPositionMap.size()];
        if(this.cacheBitSetFingerprint != null && this.cacheListToGenerateBitFingerprint.size() == aListOfUniqueSmiles.size()) {
            Collections.sort(aListOfUniqueSmiles);
            Collections.sort(this.cacheListToGenerateBitFingerprint);
            if(aListOfUniqueSmiles.equals(this.cacheListToGenerateBitFingerprint)) {
                for (int tmpPositivePositions : this.cacheBitSetFingerprint.getSetbits()) {
                    tmpBitArray[tmpPositivePositions] = 1;
                }
                //return tmpBitArray;
            }
        } else {
            this.cacheBitSetFingerprint = null;
            this.cacheListToGenerateBitFingerprint = null;
            this.getBitFingerprint(aListOfUniqueSmiles);
            for (int tmpPositivePositions : this.cacheBitSetFingerprint.getSetbits()) {
                tmpBitArray[tmpPositivePositions] = 1;
            }
        }
        return tmpBitArray;
    }
    //

    /**
     * Generates count array (float values) for the specified list (molecule).
     * Among other things, already generated results are used to generate the array.
     * For example, if a count fingerprint has already been generated for the given list of unique SMILES or for
     * the given molecule, the result of the count fingerprint is expanded into an array. Otherwise,
     * the count fingerprint is generated first and then the count array.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * @return float[] count array
     */
    private float[] createCountFloatArray(List<String> aListOfUniqueSmiles) {
        float[] tmpCountFloatArray = new float[this.uniqueSmilesToPositionMap.size()];
        if (this.cacheRawCountMap != null && this.cacheListToGenerateCountFingerprint.size() == aListOfUniqueSmiles.size()) {
            Collections.sort(aListOfUniqueSmiles);
            Collections.sort(this.cacheListToGenerateCountFingerprint);
            if (aListOfUniqueSmiles.equals(this.cacheListToGenerateCountFingerprint)) {
                for (int tmpPositivePositions: this.cacheRawCountMap.keySet()) {
                    //create float variant of cacheRawCountMap
                    tmpCountFloatArray[tmpPositivePositions] = this.cacheRawCountMap.get(tmpPositivePositions);
                }
            }
        } else {
            this.cacheRawCountMap = null;
            this.cacheListToGenerateCountFingerprint = null;
            this.getCountFingerprint(aListOfUniqueSmiles);
            for (int tmpPositivePositions: this.cacheRawCountMap.keySet()) {
                //same float variant applies here for cacheRawCountMap
                tmpCountFloatArray[tmpPositivePositions] = this.cacheRawCountMap.get(tmpPositivePositions);
            }
        }
        return tmpCountFloatArray;
    }
    //

    /**
     * Generates bit array (float values) for the specified list (molecule).
     * Among other things, already generated results are used to generate the array.
     * For example, if a bit fingerprint has already been generated for the given list of unique SMILES or for
     * the given molecule, the result of the bit fingerprint is expanded into an array. Otherwise,
     * the bit fingerprint is generated first and then the bit array.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * @return float[] bit array
     */
    private float[] createFloatBitArray(List<String> aListOfUniqueSmiles, BitSetFingerprint aBitSetFingerprint) {
        float[] tmpFloatBitArray = new float[this.uniqueSmilesToPositionMap.size()];
        if(aBitSetFingerprint != null) {
            System.out.println(aBitSetFingerprint);
            //Collections.sort(aListOfUniqueSmiles);
            //Collections.sort(this.cacheListToGenerateBitFingerprint);
            //if (aListOfUniqueSmiles.equals(this.cacheListToGenerateBitFingerprint)) {
                for (int tmpPositivePositions : aBitSetFingerprint.getSetbits()) {
                    tmpFloatBitArray[tmpPositivePositions] = 1.0f;
                }
            //}
        }
//        else {
//            this.cacheBitFingerprint = null;
//            this.cacheListToGenerateBitFingerprint = null;
//            this.getBitFingerprint(aListOfUniqueSmiles);
//            for (int tmpPositivePositions: this.cacheBitFingerprint.getSetbits()) {
//                tmpFloatBitArray[tmpPositivePositions] = 1.0f;
//            }
//        }
        return tmpFloatBitArray;
    }

    /**
     * Generates a float matrix with dimensions set by the number of molecules in form of Strings given by aFragmentsUniqueSmilesList
     * and a float array with the number of fingerprints and molecule descriptors with their respective values
     * (m. all components = fingerprint components (= Bit Array) + descriptor components).
     *
     * @param aFragmentsUniqueSmilesLists is a list of lists that store fragments in the form of unique SMILES, for a number of starter molecules.
     * @param aColumnLength defines the column length of the generated matrix (how many rows there should be).
     * @param aRowLength defines the row length of the generated matrix (how many columns there should be).
     * @param anUseBitArrayStatement defines usage of bit array or count array for fingerprint generation.
     * @return float[][] matrix with fingerprints, and if wanted, more space for additional fingerprint components.
     */
    private float[][] createFragmentsComponentsFloatMatrix(
            ArrayList<List<String>> aFragmentsUniqueSmilesLists,
            int aColumnLength,
            int aRowLength,
            boolean anUseBitArrayStatement
    ) {
        //float[column size][row size]
        float[][] tmpDataMatrix = new float[aColumnLength][aRowLength];
        float[] tmpFloatFingerprint = new float[this.fragmentsForBitSetArray.length];

        //for(each molecule's fragments) {
        // create fingerprint
        // add fingerprint to matrix
        // }
        //makes no sense right now... no idea what input gets through here
        //but for-i makes sense, could essentially lead to only one iteration through matrix for generating fingerprint and putting it in matrix
        for (int i = 0; i < aColumnLength; i++) {
            //ToDo: allows possibly to parallelize matrix generation over mols (think stack.push and .pop with index where to save result (which row it is))
            //more like fingerprint generation with an overlooking matrix generation
            if (anUseBitArrayStatement) {
                tmpFloatFingerprint = this.getFloatBitFingerprint(aFragmentsUniqueSmilesLists.get(i), this.uniqueSmilesToPositionMap, tmpFloatFingerprint);
            } else {
                tmpFloatFingerprint = this.getFloatCountFingerprint(aFragmentsUniqueSmilesLists.get(i), this.uniqueSmilesToPositionMap, tmpFloatFingerprint);
            }
            tmpDataMatrix[i] = tmpFloatFingerprint;
        }
        return tmpDataMatrix;
    }

    private float[][] fillFragmentsComponentsFloatMatrix(
            ArrayList<List<String>> aFragmentsUniqueSmilesListsArrayList,
            float[][] aFloatDataMatrix,
            boolean anUseBitArrayStatement
    ) {
        //float[row count][column count]
        for (int i = 0; i < aFloatDataMatrix.length; i++) {
            if (anUseBitArrayStatement) {
                aFloatDataMatrix[i] = this.getFloatBitFingerprint(aFragmentsUniqueSmilesListsArrayList.get(i), this.uniqueSmilesToPositionMap, aFloatDataMatrix[i]);
            } else {
                aFloatDataMatrix[i] = this.getFloatCountFingerprint(aFragmentsUniqueSmilesListsArrayList.get(i), this.uniqueSmilesToPositionMap, aFloatDataMatrix[i]);
            }
        }
        return aFloatDataMatrix;
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
        for(String tmpKeyFragmentPositionInFingerprint : this.uniqueSmilesToPositionMap.keySet()) {
            tmpPredefinedFragmentsInArray[this.uniqueSmilesToPositionMap.get(tmpKeyFragmentPositionInFingerprint)] = tmpKeyFragmentPositionInFingerprint;
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
            if (tmpUniqueSmiles.isBlank() || tmpUniqueSmiles.isEmpty()) {
                throw new IllegalArgumentException(anArgumentElementBlankEmptyMessage);
            }
        }
    }
    //
    /**
     * The key fragments defined during initialization are stored in a map.
     * The map maps the key fragments to their positions in the fingerprint.
     *
     * @return HashMap<String,Integer>
     */
    @Deprecated
    private HashMap<String,Integer> buildUniqueSmilesToPositionMap() {
        this.uniqueSmilesToPositionMap = new HashMap<>((int) (this.fragmentsForBitSetArray.length * 1.5f), 0.75f);
        int tmpValuePosition = 0;
        for (String tmpKey : this.fragmentsForBitSetArray) {
            if (!this.uniqueSmilesToPositionMap.containsKey(tmpKey)) {
                this.uniqueSmilesToPositionMap.put(tmpKey, tmpValuePosition);
                tmpValuePosition++;
            } else {
                continue;
            }
        }
        return this.uniqueSmilesToPositionMap;
    }
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
