#!/bin/bash
echo "...OK. I TIDY UP vast_out"

echo "I'm tidying up the GO files... if there's any"
mkdir vast_out/GO_ENGIDs_lists
mv vast_out/IR_DOWN-*.txt vast_out/GO_ENGIDs_lists/
mv vast_out/IR_UP-*.txt vast_out/GO_ENGIDs_lists/
mv vast_out/BG-*.txt vast_out/GO_ENGIDs_lists/
mv vast_out/AltEx-*.txt vast_out/GO_ENGIDs_lists/
mv vast_out/All_Ev-*.txt vast_out/GO_ENGIDs_lists/
rmdir vast_out/GO_ENGIDs_lists

echo "I'm tidying up the DiffAS files with dPSI... if there's any"
mkdir vast_out/DiffAS_dPSI0_files_with_dPSI
mv vast_out/CR-*-dPSI0-range*-with_dPSI.tab vast_out/DiffAS_dPSI00_files_with_dPSI/
mv vast_out/CS-*-dPSI0-range*-with_dPSI.tab vast_out/DiffAS_dPSI00_files_with_dPSI/
mv vast_out/AS_NC-*-dPSI0-range*-dPSI10-*-with_dPSI*.tab vast_out/DiffAS_dPSI0_files_with_dPSI/
mv vast_out/DiffAS-*-dPSI0-range*-with_dPSI.tab vast_out/DiffAS_dPSI0_files_with_dPSI/
rmdir vast_out/DiffAS_dPSI0_files_with_dPSI
mkdir vast_out/DiffAS_dPSI10_files_with_dPSI
mv vast_out/CR-*-dPSI10-range*-with_dPSI.tab vast_out/DiffAS_dPSI10_files_with_dPSI/
mv vast_out/CS-*-dPSI10-range*-with_dPSI.tab vast_out/DiffAS_dPSI10_files_with_dPSI/
mv vast_out/AS_NC-*-dPSI10-range*-dPSI10-*-with_dPSI*.tab vast_out/DiffAS_dPSI10_files_with_dPSI/
mv vast_out/DiffAS-*-dPSI10-range*-with_dPSI.tab vast_out/DiffAS_dPSI10_files_with_dPSI/
rmdir vast_out/DiffAS_dPSI10_files_with_dPSI
mkdir vast_out/DiffAS_dPSI15_files_with_dPSI
mv vast_out/CR-*-dPSI15-range*-with_dPSI.tab vast_out/DiffAS_dPSI15_files_with_dPSI/
mv vast_out/CS-*-dPSI15-range*-with_dPSI.tab vast_out/DiffAS_dPSI15_files_with_dPSI/
mv vast_out/AS_NC-*-dPSI15-range*-with_dPSI*.tab vast_out/DiffAS_dPSI15_files_with_dPSI/
mv vast_out/DiffAS-*-dPSI15-range*-with_dPSI.tab vast_out/DiffAS_dPSI15_files_with_dPSI/
rmdir vast_out/DiffAS_dPSI15_files_with_dPSI

echo "I'm tidying up the other DiffAS files... if there's any"
mkdir vast_out/DiffAS_dPSI10_files
mv vast_out/CR-*-dPSI10-range*.tab vast_out/DiffAS_dPSI10_files/
mv vast_out/CS-*-dPSI10-range*.tab vast_out/DiffAS_dPSI10_files/
mv vast_out/AS_NC-*-dPSI10-range*.tab vast_out/DiffAS_dPSI10_files/
mv vast_out/DiffAS-*-dPSI10-range*.tab vast_out/DiffAS_dPSI10_files/
rmdir vast_out/DiffAS_dPSI10_files
mkdir vast_out/DiffAS_dPSI0_files
mv vast_out/CR-*-dPSI0-range*.tab vast_out/DiffAS_dPSI0_files/
mv vast_out/CS-*-dPSI0-range*.tab vast_out/DiffAS_dPSI0_files/
mv vast_out/AS_NC-*-dPSI0-range*.tab vast_out/DiffAS_dPSI0_files/
mv vast_out/DiffAS-*-dPSI0-range*.tab vast_out/DiffAS_dPSI0_files/
rmdir vast_out/DiffAS_dPSI0_files
mkdir vast_out/DiffAS_dPSI15_files
mv vast_out/CR-*-dPSI15-range*.tab vast_out/DiffAS_dPSI15_files/
mv vast_out/CS-*-dPSI15-range*.tab vast_out/DiffAS_dPSI15_files/
mv vast_out/AS_NC-*-dPSI15-range*.tab vast_out/DiffAS_dPSI15_files/
mv vast_out/DiffAS-*-dPSI15-range*.tab vast_out/DiffAS_dPSI15_files/
rmdir vast_out/DiffAS_dPSI15_files
echo "I'm done. Enjoy!"
