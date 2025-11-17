$lines = Get-Content 'ProgramPhasing.cpp'
$start = 0
$end = 0
for ($i = 0; $i -lt $lines.Count; $i++) {
    if ($lines[$i] -match '^int predict_phasing_without_parents\(') {
        $start = $i
    }
    if ($start -gt 0 -and $lines[$i] -match '^int predict_phasing_with_parents_providing_GT\(') {
        $end = $i
        break
    }
}
$func1 = $lines[$start..($end-1)] -join "`n"
$func1 = $func1 -replace 'IDjob', '*IDjob'
$func1 = $func1 -replace 'placefirttoconsider', '*placefirttoconsider'
$func1 | Out-File -Encoding utf8 'temp_func1.txt'

$start2 = 0
$end2 = 0
for ($i = 0; $i -lt $lines.Count; $i++) {
    if ($lines[$i] -match '^int predict_phasing_with_parents_providing_GT\(') {
        $start2 = $i
    }
    if ($start2 -gt 0 -and $i -gt $start2 -and $lines[$i] -match '^[a-zA-Z_].*\(.*\)\s*\{' -and $lines[$i] -notmatch 'predict_phasing') {
        $end2 = $i
        break
    }
}
if ($end2 -eq 0) { $end2 = $lines.Count - 1 }
$func2 = $lines[$start2..$end2] -join "`n"
$func2 = $func2 -replace 'IDjob', '*IDjob'
$func2 = $func2 -replace 'placefirttoconsider', '*placefirttoconsider'
$func2 | Out-File -Encoding utf8 'temp_func2.txt'

Write-Host "Extraction complete. func1: $($start+1) to $end, func2: $($start2+1) to $end2"


