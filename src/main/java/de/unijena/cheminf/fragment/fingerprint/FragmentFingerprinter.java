/*
 * MIT License
 *
 * Copyright (c) 2022 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
import org.openscience.cdk.exception.UnsupportedChemObjectException;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Class to generate fingerprints. The fingerprints that can be generated are the bit fingerprint and the count fingerprint
 * Fragment fingerprints are  key-based fingerprints. Thus, the class requires to be predefined structures/fragments in
 * the form of unique SMILES to create the fingerprint. These structures must be passed when the
 * class is called (in the constructor). The class implements the interface IFragmentFingerprinter,
 * which inherits the IFingerprinter (CDK), which allows the class to compute fingerprints in 2 ways.
 * The first way to calculate a bit or count fingerprint is to perform a substructure comparison with all
 * predefined fragments for a given IAtomContainer. The second way to calculate fingerprints is by comparing
 * given fragments, which are in the form of "unique SMILES", with the predefined fragments.
 * The second possibility is thus based on a pure comparison of strings.
 *
 * @author Betuel Sevindik
 */
public class FragmentFingerprinter implements IFragmentFingerprinter {
    //<editor-fold desc="private  class variables" defaultstate="collapsed">
    /**
     * Is the array containing all the unique predefined SMILES fragments based
     * on which the fingerprints are then created.
     *
     */
    String[] fragmentArray;
    /**
     * The fragmentArray is converted into a HashMap to speed up the matching of the unique SMILES.
     * The Map maps the unique SMILES of the predefined fragments to the position they have in the input list.
     *
     */
    HashMap<String, Integer> uniqueSmilesToPositionMap;
    //</editor-fold>
    //

    /**
     * Constructor.
     * Initialization of the fingerprinter by using a user-defined set of fragments in the form of unique SMILES.
     * This list of predefined fragments is converted into a HashMap that maps the predefined fragments
     * to the position the fragments have in the list aFragmentList.
     *
     * @param aFragmentList is the ist in which the predefined fragments are stored.
     * @throws NullPointerException is thrown if the list aFragmentList is null or contains blank strings.
     */
    public FragmentFingerprinter(List<String> aFragmentList) throws NullPointerException {
        // Check whether aFragmentList is null or whether there are elements (strings) in the list that are empty.
        Objects.requireNonNull(aFragmentList, "aFragmentList (list of string instances) is null.");
        for(String tmpDefinedUniqueSMILES : aFragmentList) {
            Objects.requireNonNull(tmpDefinedUniqueSMILES, "aFragmentList (at least one list element) is null.");
            if(tmpDefinedUniqueSMILES.isBlank()) {
                throw new NullPointerException("aFragmentList (at least one list element) is blank.");
            }
        }
        this.fragmentArray = aFragmentList.toArray(new String[0]);
        this.uniqueSmilesToPositionMap = new HashMap<>(this.fragmentArray.length);
        int tmpValuePosition = 0;
        for (String tmpKey : this.fragmentArray) {
            this.uniqueSmilesToPositionMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
    }
    public FragmentFingerprinter() {

    }
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">

    /**
     * Method to generate the bit fingerprint.
     * An entered list of unique SMILES is compared with the predefined fragments.
     * If there is a match, the position of the unique SMILES is determined from the map and set to true in the
     * initialized BitSet.
     *
     * @param aListOfUniqueSmiles is a list that stores fragments in the form of unique SMILES.
     * To be able to calculate the fingerprint for a molecule, the fragments should belong to one molecule.
     * @return BitSet. BitSet is a CDK class that implements the IBitFingerprint interface of CDK.
     * This allows methods to be used that return useful information from the calculated bit fingerprint, such as the number of positive bits in the fingerprint, etc.
     * @throws NullPointerException is thrown if the list aListOfuniqueSmiles is null or contains blank strings.
     */
    @Override
    public IBitFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles) throws NullPointerException {
        Objects.requireNonNull(aListOfUniqueSmiles, "aFragmentList (list of string instances) is null.");
        for(String tmpSmiles : aListOfUniqueSmiles) {
            Objects.requireNonNull(tmpSmiles, "aFragmentList (at least one list element) is null.");
            if(tmpSmiles.isBlank()) {
                throw new NullPointerException("aFragmentList (at least one list element) is blank.");
            }
        }
        BitSet tmpBitSet = new BitSet(this.fragmentArray.length);
        for (String tmpUniqueSmiles : aListOfUniqueSmiles) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = uniqueSmilesToPositionMap.get(tmpUniqueSmiles);
                tmpBitSet.set(tmpPosition,true);
            }
        }
        return new BitSetFingerprint(tmpBitSet);
    }

    /**
     * Method to generate count fingerprint.
     * An input map that assigns unique SMILES to, the frequency of these unique SMILES is compared with the predefined
     * fragments. In case of a match, the corresponding position of the predefined fragment or the unique SMILES
     * is determined and additionally the frequency of these unique SMILES is also determined from the input map
     * and stored in a locally initialized map. I.e. this HashMap maps the position of the unique SMILES to the
     * frequency of this unique SMILES, which is specified with the input map.
     *
     * @param aUniqueSmilesToFrequencyMap  is a map that maps fragments in the form of unique SMILES to the frequency of unique SMILES.
     * To be able to calculate the fingerprint for a molecule, the fragments must belong to a molecule.
     * @return count fingerprint
     * @throws NullPointerException  is thrown if the map aUniqueSmilesToFrequencyMap is null or contains keys or values that are null or blank, respectively.
     */
    @Override
    public ICountFingerprint getCountFingerprint(Map<String, Integer> aUniqueSmilesToFrequencyMap) throws NullPointerException {
        HashMap<Integer, Integer> tmpRawMap = new HashMap<>(this.fragmentArray.length);
        Objects.requireNonNull(aUniqueSmilesToFrequencyMap, "aUniqueSmilesToFrequencyMap (Map of string and integer instances) is null.");
        for(String tmpKeys :aUniqueSmilesToFrequencyMap.keySet()) {
            if(tmpKeys == null || aUniqueSmilesToFrequencyMap.get(tmpKeys) == null || tmpKeys.isBlank()) {
                throw new NullPointerException("aUniqueSmilesToFrequencyMap (Map of string and integer instances) contains " +
                        "instances that are null or strings that are blank.");
            }
        }
        for (String tmpUniqueSmiles : aUniqueSmilesToFrequencyMap.keySet()) {
            if (this.uniqueSmilesToPositionMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = this.uniqueSmilesToPositionMap.get(tmpUniqueSmiles);
                tmpRawMap.put(tmpPosition,aUniqueSmilesToFrequencyMap.get(tmpUniqueSmiles));
               // this.countArray[tmpPosition] = aMoleculeFragmentMap.get(tmpUniqueSmiles);
            }
        }
      return new CountFingerprint(this.fragmentArray, this.uniqueSmilesToPositionMap, tmpRawMap);
    }

    /**
     * Method to generate count fingerprint.
     * The method works the same as the getCountFingerprint() above.
     * The only difference is that a list is entered as input.
     * A unique SMILES can occur more than once in the list. This list is converted into a HashMap.
     * The HashMap maps the unique SMILES to the number of occurrences of these unique SMILES in the list.
     *
     * @param aUniqueSmilesToFrequencyList is a list that stores fragments in the form of unique SMILES.
     * If a fragment occurs more than once in the molecule, it is also present more than
     * once in the list. To be able to calculate the fingerprint for a molecule,
     * the fragments should belong to one molecule.
     * @return count fingerprint
     * @throws NullPointerException
     */
    @Override
    public ICountFingerprint getCountFingerprint(List<String>  aUniqueSmilesToFrequencyList) throws NullPointerException {
        HashMap<String, Integer> tmpUniqueSmilesToFrequencyCountMap = new HashMap<>();
        int i = 1;
        Objects.requireNonNull(aUniqueSmilesToFrequencyList, "aUniqueSmilesToFrequencyList (list of string instances) is null.");
        for (String tmpSmiles :  aUniqueSmilesToFrequencyList) {
            Objects.requireNonNull(tmpSmiles, "aUniqueSmilesToFrequencyList (at least one list element) is null.");
            if(tmpSmiles.isBlank()) {
                throw new NullPointerException("aUniqueSmilesToFrequencyList (at least one list element) is blank.");
            }
            if (tmpUniqueSmilesToFrequencyCountMap.containsKey(tmpSmiles) == false) {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, i);
            } else {
                tmpUniqueSmilesToFrequencyCountMap.put(tmpSmiles, tmpUniqueSmilesToFrequencyCountMap.get(tmpSmiles) + 1);

            }
        }
        return this.getCountFingerprint(tmpUniqueSmilesToFrequencyCountMap);
    }

    /**
     *
     * @return
     */
    @Override
    public String getVersionDescription() {
       throw new UnsupportedOperationException();
    }

    /**
     *
     * @param mol molecule
     * @return
     * @throws CDKException
     */
    @Override
    public BitSet getFingerprint(IAtomContainer mol) throws CDKException {
        throw new UnsupportedOperationException();
    }

    /**
     *
     * @param container {@link IAtomContainer} for which the fingerprint should be calculated.
     * @return
     * @throws CDKException
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

    /**
     *
     * @param container {@link IAtomContainer} for which the fingerprint should be calculated.
     * @return
     * @throws CDKException
     */
    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

    /**
     *
     * @param container IAtomContainer for which the fingerprint should be calculated.
     * @return
     * @throws CDKException
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {
       throw new UnsupportedOperationException();
    }

    /**
     * Since the FragmentFingerprinter is a key-based fingerprint, the size of the fingerprint
     * corresponds to the number of predefined fragments(unique SMILES)
     *
     * @return int
     */
    @Override
    public int getSize() {
     return this.fragmentArray.length;
    }

    /**
     *  returns the bit definitions i.e. which  bit stands for which fragment SMILES
     *
     * @param aBit
     * @return unique SMILES
     */
    public String getBitDefinition(int aBit) {
      return this.fragmentArray[aBit];
    }
    public int[] getBitArray(){
        throw new UnsupportedOperationException();
    }
    public int[] getCountArray(){
        throw new UnsupportedOperationException();
    }
}
