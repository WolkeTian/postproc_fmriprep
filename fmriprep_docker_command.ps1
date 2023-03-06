# 开始计时
$sw = [Diagnostics.Stopwatch]::StartNew()

$subids_1 = '01', '02', '03', '04', '05', '06', '07', '08', '09'
$subids_2 = [string] 10..50

$subids = $subids_1 + $subids_2

foreach ($indexs in 0..49) 
{
$proc_subs = $subids($indexs) # 索引从0开始
$proc_subs

docker run --rm -it  `
-v E:\Projects\SD_PANAS\fmriprep_datasets:/input `
-v E:\Projects\SD_PANAS\prep_out:/output `
-v D:\google_drive_local\archive\FreesurferLicense\license.txt:/license `
nipreps/fmriprep:latest /input /output participant --fs-license-file /license `
--skip-bids-validation --use-aroma --cifti-output `
--output-spaces MNI152NLin2009cAsym fsnative fsaverage --notrack --mem_mb 24000 `
--participant_label $proc_subs
}

# 停止計時
$sw.Stop()

# 輸出測量時間
$sw.Elapsed
