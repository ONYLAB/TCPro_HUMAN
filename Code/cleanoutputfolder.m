folders = readtable('folders.dat','Delimiter','\t');
for i = 1:height(folders)
drug = folders.folders{i};
cd(drug)
cd Output
unix('rm *')
cd ..
cd ..
end