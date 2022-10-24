package de.unijena.cheminf.fragment.fingerprint;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class CountFragmentFingerprint implements IFingerprint {

    public CountFragmentFingerprint() {

    }

    @Override
    public ArrayList generateFragmentFingerprint(ArrayList<String> aFragmentList, HashMap<String, String> aMoleculeFragments) {
        int tmpVectorSize = aFragmentList.size();
        int[] tmpCountVector = new int[tmpVectorSize];
        ArrayList<Object> yeni = new ArrayList<>();
        ArrayList<Integer> fgh = new ArrayList<>();
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            tmpCountVector[tmpDefaultBit] = 0;
        }
        HashMap<String, Integer> tmpFragmentHashMap = new HashMap<>();
        // sort aFragmentList alphabetically
        Collections.sort(aFragmentList, String.CASE_INSENSITIVE_ORDER);
        int tmpValuePosition = 0;
        for (String tmpKey : aFragmentList) {
            tmpFragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
        for (String str : aMoleculeFragments.keySet()) {
            if (tmpFragmentHashMap.containsKey(str) == true) {
                int tmpPosition = tmpFragmentHashMap.get(str);
                fgh.get(tmpPosition);
                tmpCountVector[tmpPosition] = Integer.parseInt(aMoleculeFragments.get(str));
            }
        }
        yeni.add(tmpCountVector);
        yeni.add(fgh);
        return yeni;
       // return tmpCountVector;


    }

    @Override
    public int getFragmentFingerprintSize() {
        return 0;
    }

    @Override
    public Object getNumberPositiveBits(ArrayList<Object> aList) {
        return aList.get(1);
    }

    @Override
    public int getIndicesPositiveBits() {
        return 0;
    }
}
