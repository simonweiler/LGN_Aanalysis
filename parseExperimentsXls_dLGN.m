function [batchopt] = parseExperimentsXls_dLGN(path,user)

[xls_num,xls_txt]=xlsread(path);

loadcol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'BatchAnalyze')));
mousecol      = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExperimentalDay')));
expcol        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsSW')));
expcols2        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsSW_slice_nr')));
if user==1
expcol      = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsMF')));
expcols2    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsMF_slice_nr')));
end
% samesitecol   = find(~cellfun(@isempty, strfind(xls_txt(1,:),'SameSiteID')));
% baselinecol   = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Baseline')));
% basepaircol   = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Baselinepair')));
% basepair14col = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Baselinepair14')));
% recoverycol   = find(~cellfun(@isempty, strfind(xls_txt(1,:),'fullRec')));
 loaddrivecol  = find(~cellfun(@isempty, strfind(xls_txt(1,:),'loaddrive')));
 animalname    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Animal_ID')));
% airpuffcol    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Airpuff')));
% spontcol      = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExpIDsSpont')));
% sftfcol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExpIDsSFTF')));

k = 1;

batchopt.XLS.txt = xls_txt;
batchopt.XLS.num = xls_txt;

for i = 2:size(xls_txt,1);
    ana{k}= xls_num(i-1,loadcol-1);
    
    
    
    if ~ana{k}
        disp(['skipping experiments ' xls_txt(i,mousecol) '(no batchload flag)']);
        continue
    end
    
    batchopt.mouse{k} = xls_txt(i,mousecol);
    batchopt.mouseID{k} = xls_txt(i,animalname);
    
    expcellids{k}                = xls_txt(i,expcol);
     expcellids2{k}                = xls_txt(i,expcols2);
%     spontcellids{k}                = xls_txt(i,spontcol);
%     sftfcellids{k}                = xls_txt(i,sftfcol);
%     puffcellids{k}                = xls_txt(i,airpuffcol);
    
    %batchopt.exp_ids{k}          = (expcellids{k}{1});


    batchopt.exp_ids{k}          = str2num((expcellids{k}{1}));
    batchopt.exp_ids2{k}          = str2num((expcellids2{k}{1}));
%     batchopt.spont_ids{k}          = str2num((spontcellids{k}{1}));
%     batchopt.sftf_ids{k}          = str2num((sftfcellids{k}{1}));
%     batchopt.puff_ids{k}          = str2num((puffcellids{k}{1}));
%     
%     batchopt.samesite{k}         = xls_num(i-1,samesitecol-1);
%     batchopt.baseline{k}         = xls_num(i-1,baselinecol(1)-1);
%     try
%         batchopt.baselinepair{k}     = eval(cell2mat(xls_txt(i,basepaircol(1))));
%         batchopt.baselinepair14{k}   = eval(cell2mat(xls_txt(i,basepair14col)));
%         batchopt.recovery{k}         = xls_num(i-1,recoverycol-1);
%     end
    batchopt.loaddrive{k}        = xls_txt(i,loaddrivecol);
    k = k+1;
end

