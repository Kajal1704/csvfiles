close all;
clear all;

srcFiles = dir('D:\HD paper - dataset\*.fasta');
W=zeros(25,17); %%%% (zeros(no. of i/p samples, no.of feature for one sample))
for m = 1 : length(srcFiles)  
    m
    
    filename = strcat('D:\HD paper - dataset\',srcFiles(m).name);
    
[header,seqs] = fastaread(filename)

%seqs='CTCCGAGACG';

len=length(seqs)


%%%%%% ASCII values representation %%%%
tic
as_val=double(seqs)' %%% ASCII values
toc
 
%%%%% Feature extraction %%%%%%%%%%
%glcm=graycomatrix(z1,'Offset',[0 1;1 0]);

glcm1 = graycomatrix(as_val,'Offset',[0 1]) %%% default offset [0 1], Angle 0 Degree
features1 = graycoprops(glcm1,{'contrast','homogeneity','energy','correlation'});
contrast1 = features1.Contrast;
homogeneity1 = features1.Homogeneity;
energy1 = features1.Energy;
correlation1 = abs(features1.Correlation);
f1_0_Deg = [contrast1,homogeneity1,energy1,correlation1];  %%% f1_0_Deg (i.e feature1 zero Degree)

glcm2 = graycomatrix(as_val,'Offset',[1 0]);  %%% offset [1 0], Angle 180 Degree
features2 = graycoprops(glcm2,{'contrast','homogeneity','energy','correlation'});
contrast2 = features2.Contrast;
homogeneity2 = features2.Homogeneity;
energy2 = features2.Energy;
correlation2 = abs(features2.Correlation);
f2_180_Deg = [contrast2,homogeneity2,energy2,correlation2];   %%% f2_180_Deg (i.e feature2 180 Degree)


glcm3 = graycomatrix(as_val,'Offset',[0 0]);  
features3 = graycoprops(glcm3,{'contrast','homogeneity','energy','correlation'});
contrast3 = features3.Contrast;
homogeneity3 = features3.Homogeneity;
energy3 = features3.Energy;
correlation3 = abs(features3.Correlation);
f3_00_No_Deg = [contrast3,homogeneity3,energy3,correlation3]; %%% f3_00_No_Deg (i.e feature3 00 No Degree)

glcm4 = graycomatrix(as_val,'Offset',[-1 1]); %%% offset [-1 1], Angle 45 Degree
features4 = graycoprops(glcm4,{'contrast','homogeneity','energy','correlation'});
contrast4 = features4.Contrast;
homogeneity4 = features4.Homogeneity;
energy4 = features4.Energy;
correlation4 = abs(features4.Correlation);
f4_45_Deg = [contrast4,homogeneity4,energy4,correlation4];  %%% f4_45_Deg (i.e feature4 45 Degree)

glcm5 = graycomatrix(as_val,'Offset',[-1 0]); %%% offset [-1 0], Angle 90 Degree
features5 = graycoprops(glcm5,{'contrast','homogeneity','energy','correlation'});
contrast5 = features5.Contrast;
homogeneity5 = features5.Homogeneity;
energy5 = features5.Energy;
correlation5 = abs(features5.Correlation);
f5_90_Deg = [contrast5,homogeneity5,energy5,correlation5];   %%% f5_90_Deg (i.e feature5 90 Degree)


glcm6 = graycomatrix(as_val,'Offset',[-1 -1]); %%% offset [-1 -1], Angle 135 Degree
features6 = graycoprops(glcm6,{'contrast','homogeneity','energy','correlation'});
contrast6 = features6.Contrast;
homogeneity6 = features6.Homogeneity;
energy6 = features6.Energy;
correlation6 = abs(features6.Correlation);
f6_135_Deg = [contrast6,homogeneity6,energy6,correlation6];  %%% f6_135_Deg (i.e feature6 135 Degree)

glcm7 = graycomatrix(as_val,'Offset',[1 -1]); %%% offset [1 -1], Angle 225 Degree
features7 = graycoprops(glcm7,{'contrast','homogeneity','energy','correlation'});
contrast7 = features7.Contrast;
homogeneity7 = features7.Homogeneity;
energy7 = features7.Energy;
correlation7 = abs(features7.Correlation);
f7_225_Deg = [contrast7,homogeneity7,energy7,correlation7];  %%% f7_225_Deg (i.e feature7 225 Degree)

glcm8 = graycomatrix(as_val,'Offset',[0 -1]); %%% offset [0 -1], Angle 270 Degree
features8 = graycoprops(glcm8,{'contrast','homogeneity','energy','correlation'});
contrast8 = features8.Contrast;
homogeneity8 = features8.Homogeneity;
energy8 = features8.Energy;
correlation8 = abs(features8.Correlation);
f8_270_Deg = [contrast8,homogeneity8,energy8,correlation8];  %%% f8_270_Deg (i.e feature8 270 Degree)

glcm9 = graycomatrix(as_val,'Offset',[1 1]); %%% offset [1 1], Angle 315 Degree
features9 = graycoprops(glcm9,{'contrast','homogeneity','energy','correlation'});
contrast9 = features9.Contrast;
homogeneity9 = features9.Homogeneity;
energy9 = features9.Energy;
correlation9 = abs(features9.Correlation);
f9_315_Deg = [contrast9,homogeneity9,energy9,correlation9];  %%% f9_315_Deg (i.e feature9 315 Degree)

GLCM_features=[f1_0_Deg' f2_180_Deg' f3_00_No_Deg' f4_45_Deg' f5_90_Deg' f6_135_Deg' f7_225_Deg' f8_270_Deg' f9_315_Deg'];

z1_ent=entropy(as_val);
z1_mean=mean(as_val);
z1_median=median(as_val);
z1_stand_dev=std(as_val);
bases_cnt=basecount(seqs);

%%%% CAG repeat count %%%%
l=length(seqs);
kmer=nmercount(seqs,3)

kmer_strings=kmer(:,1);
kmer_counts=kmer(:,2);

l_kmer_strings=length(kmer_strings);
l_kmer_counts=length(kmer_counts);

CAG_repeat=0;
for i=1:l_kmer_strings
    if(kmer_strings{i}=='CAG')
         CAG_repeat=CAG_repeat+cell2mat(kmer_counts(i));
         %disp(CAG_repeat)
     end
     %elseif
         %CAG_notrepeat=0;
         %disp(CAG_notrepeat)
         %CAG_repeat={0};
         %disp(CAG_notrepeat)
     end
%end


%CAG_count=cell2mat(CAG_repeat)   %%% cell2mat display the value as in decimal  e.g. 10.0000 
CAG_count=CAG_repeat   %%%%% this will display the value as in cell e.g. [10]


%W(i,1:20)=[struct2cell(features2) struct2cell(features3) struct2cell(features5) z1_ent z1_mean z1_median z1_stand_dev bases_cnt.A bases_cnt.C bases_cnt.G bases_cnt.T]
m
W(m,1:17)=[contrast2 homogeneity2 energy2 correlation2 contrast3 homogeneity3 energy3 correlation3 z1_ent z1_mean z1_median z1_stand_dev bases_cnt.A bases_cnt.C bases_cnt.G bases_cnt.T CAG_count]
size_W=size(W)
end
display(W);