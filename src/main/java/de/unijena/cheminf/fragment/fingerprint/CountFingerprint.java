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

import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;

/**
 * The CountFingerprint class implements the CDK interface ICountFingerprint.
 * ICountFingerprint provides useful methods to obtain information about the calculated count fingerprint.
 * Instead of using the ICountFingerprint interface that implements the CDK class IntArrayCountFingerprint,
 * a new CountFingerprint class has been created here that also implements the ICountFingerprint interface.
 * The IntArrayCountFingerprint class assumes hashed count fingerprints, while here they are key-based count
 * fingerprints, so it is necessary to create the CountFingerprint class to treat the fingerprints as
 * key-based fingerprints. This class does not allow the creation of a counting fingerprint. It is only
 * responsible for providing methods for an already created count fingerprint to get information from the fingerprint.
 *
 * @author Betuel Sevindik
 * @version 1.0.0.0
 */
public class CountFingerprint implements ICountFingerprint {
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * The HashMap maps the predefined (key) fragments/unique SMILES to the position they have in the fragment array
     * that stores all key fragments. Thus the value of the map is the position of the
     * predefined fragments in the fingerprint. // TODO ask
     */
    private final HashMap<String,Integer> uniqueSmilesToPositionMap;
    /**
     * Value for determining the initial capacity of cards
     */
    private final int INITIAL_CAPACITY_VALUE = Math.round((4/3) + 1);
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    /**
     * changes the behavior of the fingerprint. If behaveAsBitFingerprint == true, the count fingerprint
     * behaves the same as a bit fingerprint.
     */
    private boolean behaveAsBitFingerprint;
    /**
     * change the behavior of the fingerprint. If behaveAsMergeFingerprint == true, the fingerprint behaves like
     * a merged fingerprint
     */
    private boolean behaveAsMergeFingerprint;
    private HashMap<Integer,Integer> rawMapForMergedFingerprint;
    /**
     * Result Map of a count fingerprint. When the count fingerprint is created, if there is a match between
     * the key fragments/unique SMILES and the passed fragments with the corresponding frequencies, usually
     * belonging to a molecule, the position of the key fragment assigned to it during initialization is mapped
     * to the frequency of this fragment in the molecule. So, in this map, only the positions and frequencies of
     * the fragments are stored if a match has occurred. In the following, the map can also be referred to as a raw map.
     * // TODO ask
     */
    private HashMap<Integer,Integer> uniqueSmilesPositionToFrequencyCountRawMap;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     * Initialization of the CountFingerprint.
     * The specified parameters are checked for validity.
     * Duplicate fragment SMILES (key fragments) in the given array will be ignored and not be a part of the
     * count fingerprint multiple times.
     *
     * @param anArrayOfFragments is a string array that stores all predefined (key) fragments,
     *                           which are in the form of unique SMILES. The generation of the fingerprints
     *                           is based on these fragments.
     * @param aPositionToFrequencyMap This map is a raw map and must map the positions of key fragments in the
     *                                fingerprint to their frequencies in the molecule or any set of fragments.
     * @throws NullPointerException is thrown if the arguments are null.
     * @throws IllegalArgumentException is thrown if the list anArrayOfFragments contains blank/empty strings.
     */
    public CountFingerprint(String[] anArrayOfFragments, HashMap<Integer, Integer> aPositionToFrequencyMap) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(anArrayOfFragments, "anArrayOfFragments (array of string instances) is null.");
        Objects.requireNonNull(aPositionToFrequencyMap, "aPositionToFrequencyMap is null.");
        for(String tmpKeyUniqueSMILES : anArrayOfFragments) {
            Objects.requireNonNull(tmpKeyUniqueSMILES, "anArrayOfFragments (at least one list element) is null.");
            if(tmpKeyUniqueSMILES.isBlank() || tmpKeyUniqueSMILES.isEmpty()) {
                throw new IllegalArgumentException("anArrayOfFragments (at least one list element) is blank/empty.");
            }
        }
        for (Integer tmpUniqueSmiles : aPositionToFrequencyMap.keySet()) {
            if (tmpUniqueSmiles == null || aPositionToFrequencyMap.get(tmpUniqueSmiles) == null) {
                throw new NullPointerException("aPositionToFrequencyMap (Map of integer instances) contains " +
                        "instances that are null.");
            }
        }
        String[] tmpPredefinedFragments = anArrayOfFragments;
        this.uniqueSmilesToPositionMap = new HashMap<>(tmpPredefinedFragments.length*INITIAL_CAPACITY_VALUE, 0.75f);
        int tmpValuePosition = 0;
        for (String tmpKey : tmpPredefinedFragments) {
            if(!this.uniqueSmilesToPositionMap.containsKey(tmpKey)) {
                this.uniqueSmilesToPositionMap.put(tmpKey, tmpValuePosition);
                tmpValuePosition++;
            } else {
                continue;
            }
            //  this.uniqueSmilesToPositionMap.putIfAbsent(tmpKey, tmpValuePosition);
            // tmpValuePosition++;
        }
        this.uniqueSmilesPositionToFrequencyCountRawMap = aPositionToFrequencyMap;
        this.behaveAsBitFingerprint = false;
        this.behaveAsMergeFingerprint = false;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Overridden public methods " defaultstate="collapsed">
    /**
     * {@inheritDoc}
     *
     * Since fragment fingerprints are key-based, the number of bits in the fingerprint
     * is equal to the number of predefined fragments (unique SMILES) if the list of key fragments passed during
     * initialization does not contain duplicates, otherwise the size of the fingerprint may be smaller
     * than the number of fragments passed since duplicates are removed/ignored.
     *
     */
    @Override
    public long size() {
        return this.uniqueSmilesToPositionMap.size();
    }
    //
    /**
     * {@inheritDoc}
     *
     * Fragment fingerprints are key-based fingerprints,
     * therefore the number of populated bins corresponds to the number of predefined fragments (unique SMILES).
     * If the list of key fragments passed during initialization does not contain duplicates, otherwise the size
     * of the fingerprint may be smaller than the number of fragments passed since duplicates are removed/ignored.
     *
     */
    @Override
    public int numOfPopulatedbins() {
        return this.uniqueSmilesToPositionMap.size();
    }
    //
    /**
     * {@inheritDoc}
     *
     * For an index specified within the fingerprint size, the corresponding frequency is returned.
     * If the given index is greater than the size of the fingerprint or a negative value,
     * an IllegalArgumentException is thrown.
     * The size of the fingerprint corresponds to the number of predefined fragments passed during initialization.
     * However, if duplicates are included, they are ignored and are not a part of the fingerprint multiple times.
     * The size of the fingerprint in this case is then different from the number of predefined fragments.
     *
     * @throws IllegalArgumentException is thrown if the given index does not exist in the fingerprint.
     */
    @Override
    public int getCount(int index) throws IllegalArgumentException {
        if (index >= this.uniqueSmilesToPositionMap.size() || index < 0) {
            throw new IllegalArgumentException("This position does not exist in the fingerprint (undefined state).");
        } else if (index >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(index)) {
            if (this.behaveAsBitFingerprint) {
                return 1;
            } else {
                return this.uniqueSmilesPositionToFrequencyCountRawMap.get(index);
            }
        } else {
            return 0;
        }
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since this is a key-based fingerprint, the hash value is  the position of the bin in the
     * fingerprint (no hash value is calculated).
     *
     */
    @Override
    public int getHash(int index) {
        String[] tmpPredefinedFragmentsInArray = new String[this.uniqueSmilesToPositionMap.size()];
        System.out.println(tmpPredefinedFragmentsInArray.length + "---length");
        for(String tmpKeyFragmentPositionInFingerprint : this.uniqueSmilesToPositionMap.keySet()) {
            tmpPredefinedFragmentsInArray[this.uniqueSmilesToPositionMap.get(tmpKeyFragmentPositionInFingerprint)] = tmpKeyFragmentPositionInFingerprint;
        }
        if(index >= tmpPredefinedFragmentsInArray.length || index < 0 ) {
            throw new IllegalArgumentException("This hash value/position does not exist in the fingerprint (undefined state).");
        }
        return this.uniqueSmilesToPositionMap.get(tmpPredefinedFragmentsInArray[index]);
    }
    //
    /**
     * {@inheritDoc}
     *
     * Method for merging the given fingerprint fp into a current fingerprint.
     * Merging is intended only for count fingerprints generated from the same fragment set.
     */
    @Override
    public void merge(ICountFingerprint fp) {
        this.behaveAsMergeFingerprint = true;
        this.rawMapForMergedFingerprint =  this.mergeFingerprint((CountFingerprint) fp);
    }
    //
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {
        this.behaveAsBitFingerprint = behaveAsBitFingerprint;
    }
    //
    /**
     * {@inheritDoc}
     *
     * The parameter hash is not a calculated hash value, but also corresponds to the
     * position of the bin in the fingerprint.
     *
     * @throws IllegalArgumentException is thrown if the given hash value is negative.
     */
    @Override
    public boolean hasHash(int hash) throws IllegalArgumentException {
        if(hash < this.uniqueSmilesToPositionMap.size() && hash >= 0) {
            return true;
        } else if (hash < 0) {
            throw new IllegalArgumentException("Negative values are not allowed.");
        } else {
            return false;
        }
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since the fragment fingerprint is a key-based fingerprint and the hash value therefore indicates the
     * position of the respective bin in the fingerprint, each hash value occurs only once in the fingerprint.
     * Therefore, the method returns the frequency of the fragment in the specified bin.
     *
     * @see #getCount(int)
     *
     * @throws IllegalArgumentException is thrown if the given hash value do not exist in the fingerprint.
     *
     */
    @Override
    public int getCountForHash(int hash) throws IllegalArgumentException {
        System.out.println(this.rawMapForMergedFingerprint + "----gemergtes raw map tes");
        if (hash >= this.uniqueSmilesToPositionMap.size() || hash < 0) {
            throw new IllegalArgumentException("This position does not exist in the fingerprint (undefined state).");
        } else if (hash >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(hash)) {
            if (this.behaveAsBitFingerprint) {
                return 1;
            } else {
                if(!this.behaveAsMergeFingerprint)
                return this.uniqueSmilesPositionToFrequencyCountRawMap.get(hash);
                else {
                    return this.rawMapForMergedFingerprint.get(hash);
                }
            }
        } else {
            return 0;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public method" defaultstate="collapsed">
    /**
     * Method returns the frequency in the count fingerprint for a given SMILES string.
     *
     * @param aSmiles SMILES String for which the frequency is to be returned. These SMILES should represent a key fragment.
     * @return int count value
     * @throws IllegalArgumentException is thrown if the given SMILES string is not a key fragment.
     */
    public int count(String aSmiles) throws IllegalArgumentException {
        if(!this.uniqueSmilesToPositionMap.containsKey(aSmiles)) {
            throw new IllegalArgumentException("The given SMILES string is not available");
        }
       int tmpPosition =  this.uniqueSmilesToPositionMap.get(aSmiles);
        if(this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(tmpPosition)) {
            if(this.behaveAsBitFingerprint) {
                //return this.uniqueSmilesPositionToFrequencyCountRawMap.get(tmpPosition);
                return 1;
            } else {
                return this.uniqueSmilesPositionToFrequencyCountRawMap.get(tmpPosition);
            }
        } else {
            return 0;
        }
    }
    private HashMap<Integer,Integer> mergeFingerprint(CountFingerprint aCountFingerprint) {
        if(this.uniqueSmilesToPositionMap.size() != aCountFingerprint.size()) {
            throw new IllegalArgumentException("The two fingerprints are not the same size. Is only possible with" +
                    "fingerprints that come from the same fragment set.");
        }
        HashMap<Integer, Integer> tmpRawMap = new HashMap<>(this.uniqueSmilesPositionToFrequencyCountRawMap.size()+aCountFingerprint.getRawMap().size()* INITIAL_CAPACITY_VALUE, 0.75f);
        System.out.println(this.uniqueSmilesPositionToFrequencyCountRawMap.size() +"---size");
        System.out.println(this.uniqueSmilesPositionToFrequencyCountRawMap + "---ursprungliche map");
        System.out.println(aCountFingerprint.getRawMap() + "---raw map 2");
        for( int tmpKey : aCountFingerprint.getRawMap().keySet()) {
            System.out.println(tmpKey + "---tmpkey");
            // map with frequencies = 0, so in this case it is not a raw map.
            if (this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(tmpKey)) {
                System.out.println("hallo");
                tmpRawMap.put(tmpKey, this.uniqueSmilesPositionToFrequencyCountRawMap.get(tmpKey) + aCountFingerprint.getRawMap().get(tmpKey));
            } else {
                System.out.println("nooo");
                tmpRawMap.put(tmpKey, this.uniqueSmilesPositionToFrequencyCountRawMap.get(tmpKey));
            }
        }
        return tmpRawMap;
    }
    public HashMap<Integer, Integer> getRawMap() {
        return this.uniqueSmilesPositionToFrequencyCountRawMap;
    }
    //</editor-fold>
    //
}
