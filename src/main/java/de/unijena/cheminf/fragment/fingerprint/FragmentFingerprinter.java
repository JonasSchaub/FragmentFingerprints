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
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class to generate a bit fragment fingerprint
 *
 * @author Betuel Sevindik
 */
public class FragmentFingerprinter implements IFragmentFingerprinter {
    //<editor-fold desc="private  class variables" defaultstate="collapsed">
    /**
     * bitSet
     */
    BitSet bitSet;
    /**
     * Bit array to  display the complete fingerprint
     */
    int[] bitArray;
    /**
     * Count array to display the complete fingerprint
     */
    int[] countArray;
    /**
     * Is the array that contains all the unique SMILES for the fragments, on the basis of
     * which the fingerprints are then created.
     */
    String[] fragmentArray;
    /**
     * The fragmentArray is converted into a HashMap to speed up the matching of the unique SMILES. The keys in the
     * map correspond to the unique SMILES and the values to the hash values.
     */
    HashMap<String, Integer> fragmentHashMap; // TODO rename
    /**
     * Copy of bitArray
     */
    int[] copyBitArray;
    /**
     * Copy of countArray
     */
    int[] copyCountArray;
    //</editor-fold>
    //
    /**
     * Constructor
     * Initialisation of the fingerprint
     */
    public FragmentFingerprinter(List<String> aFragmentList) {
        this.fragmentArray = aFragmentList.toArray(new String[0]);
        int tmpVectorSize = this.fragmentArray.length;
        this.bitArray = new int[tmpVectorSize];
        this.copyBitArray = this.bitArray;
        this.countArray = new int[tmpVectorSize];
        this.copyCountArray = this.countArray;
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            this.bitArray[tmpDefaultBit] = 0;
            this.countArray[tmpDefaultBit] = 0;
        }
        this.fragmentHashMap = new HashMap<>();
        // sort aFragmentList alphabetically
       // Collections.sort(this.fragmentSetList, String.CASE_INSENSITIVE_ORDER);
        Arrays.sort(this.fragmentArray, String.CASE_INSENSITIVE_ORDER);
        int tmpValuePosition = 0;
        for (String tmpKey : this.fragmentArray) {
            fragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
    }
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Method to generate the bit fingerprint.
     * The fingerprint is only created by matching the unique SMILES (String matching).
     *
     * @param aListOfUniqueSmiles a list containing all the fragments produced by the fragmentation of a given molecule
     * @return IBitFingerprint (cdk)
     */
    @Override
    public IBitFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles) {
        this.bitSet = new BitSet(this.fragmentArray.length);
        for (String tmpUniqueSmiles : aListOfUniqueSmiles) {
            if (this.fragmentHashMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = fragmentHashMap.get(tmpUniqueSmiles);
                this.bitArray[tmpPosition] = 1;
                this.bitSet.set(tmpPosition,true);
            }
        }
        return new BitSetFingerprint(this.bitSet);
    }

    /**
     * Method to generate count fingerprint with a Map as input.
     * The fingerprint is only created by matching the unique SMILES (String matching).
     *
     * @param aMoleculeFragmentMap a HashMap containing all fragments of the molecule and the corresponding frequencies of the generated fragments.
     *   Key = unique SMILES, Value = frequency
     * @return CountFingerprint; new class to organise methods related to the count Fingerprint.
     */
    @Override
    public ICountFingerprint getCountFingerprint(Map<String, Integer> aMoleculeFragmentMap) {
        for (String tmpUniqueSmiles : aMoleculeFragmentMap.keySet()) {
            if (this.fragmentHashMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = this.fragmentHashMap.get(tmpUniqueSmiles);
                this.countArray[tmpPosition] = aMoleculeFragmentMap.get(tmpUniqueSmiles);
            }
        }
      return new CountFingerprint(this.fragmentArray, this.countArray, this.fragmentHashMap);
    }

    /**
     * Method to generate count fingerprint with a List as input.
     * If a unique SMILES appears more than once in the list, its frequency is summarised.
     *
     * @param aListToCountMoleculeFragmentSmiles
     * @return ICountFingerprint (cdk)
     * @author TODO add author
     *
     */
    @Override
    public ICountFingerprint getCountFingerprint(List<String> aListToCountMoleculeFragmentSmiles) {
        HashMap<String, Integer> tmpCountMap = new HashMap<>();
        int tmpUniqueSmilesFrequencyInList = 0;
        String tmpUniqueSmilesInList = null;
        for (String tmpSmiles : aListToCountMoleculeFragmentSmiles) {
            if (tmpUniqueSmilesInList == null || !tmpSmiles.equals(tmpUniqueSmilesInList)) {
                if (tmpUniqueSmilesFrequencyInList > 0)
                    tmpCountMap.put(tmpUniqueSmilesInList, tmpUniqueSmilesFrequencyInList);
                tmpUniqueSmilesFrequencyInList = 1;
                tmpUniqueSmilesInList = tmpSmiles;
            } else {
                ++tmpUniqueSmilesFrequencyInList;
            }
        }
        if (tmpUniqueSmilesFrequencyInList > 0)
            tmpCountMap.put(tmpUniqueSmilesInList, tmpUniqueSmilesFrequencyInList);
         return this.getCountFingerprint(tmpCountMap);
       // return new CountFingerprint(this.fragmentArray, this.countArray, this.fragmentHashMap);
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
     * Returns the size of the fingerprint
     *
     * @return int
     */
    @Override
    public int getSize() {
     return this.fragmentArray.length;
    }

    /**
     * Returns the complete count fingerprint
     * Input must be a map
     *
     * @return
     */
    public int[] getCountArray(Map<String, Integer> aMoleculeFragmentMap) {
        if(this.countArray == this.copyCountArray)
            this.getCountFingerprint(aMoleculeFragmentMap);
       return this.countArray;
    }

    /**
     * Returns the complete count fingerprint, but the input must be a list
     *
     * @param aListToCountMoleculeFragmentSmiles
     * @return
     */
    public int[] getCountArray(List<String> aListToCountMoleculeFragmentSmiles) {
        if(this.countArray == this.copyCountArray) {
            this.getCountFingerprint(aListToCountMoleculeFragmentSmiles);
        }
        return this.countArray;
    }

    /**
     * Returns the complete bit fingerprint
     *
     * @return
     * */
    public int[] getBitArray(List<String> aListOfUniqueSmiles) {
        if(this.bitArray == this.copyBitArray)
            this.getBitFingerprint(aListOfUniqueSmiles);
        return this.bitArray;
    }
}
