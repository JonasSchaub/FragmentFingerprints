package de.unijena.cheminf.fragment.fingerprint;



import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class BitFragmentFingerprint implements IFingerprint {
    public BitFragmentFingerprint() {

    }

    @Override
    public ArrayList generateFragmentFingerprint(ArrayList<String> aFragmentList, HashMap<String, String> aMoleculeFragments) {
        int tmpVectorSize = aFragmentList.size();
        int[] tmpBitVector = new int[tmpVectorSize];
        HashMap<int[], ArrayList<Integer>> overAllMap = new HashMap<>();
        ArrayList<Integer> fgh = new ArrayList<>();
        ArrayList<Object> yeni = new ArrayList<>();
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            tmpBitVector[tmpDefaultBit] = 0;
        }
        HashMap<String, Integer> tmpFragmentHashMap = new HashMap<>();
        // sort aFragmentList alphabetically
        Collections.sort(aFragmentList, String.CASE_INSENSITIVE_ORDER);
        System.out.println(aFragmentList + "---- sortierte fragment list");
        System.out.println(aFragmentList.size());
        int tmpValuePosition = 0;
        for (String tmpKey : aFragmentList) {
            tmpFragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;

        }
        for (String str : aMoleculeFragments.keySet()) {
            if (tmpFragmentHashMap.containsKey(str) == true) {
                int tmpPosition = tmpFragmentHashMap.get(str);
                fgh.add(tmpPosition);
                tmpBitVector[tmpPosition] = 1;

            }
        }
      //  return tmpBitVector;
        overAllMap.put(tmpBitVector,fgh);
       // return overAllMap;
        yeni.add(tmpBitVector);
        yeni.add(fgh);
        return yeni;

    }

    @Override
    public int getFragmentFingerprintSize() {
        return 0;
    }

    public Object getNumberPositiveBits(ArrayList<Object> aList) {
       return aList.get(1);

    }

    @Override
    public int getIndicesPositiveBits() {
        return 0;
    }



}
