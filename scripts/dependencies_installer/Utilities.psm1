Function PrintDecoratedMessage([string] $Name, [string] $Message) {
    $Str = "| " + $Name + " Installer: " + $Message + " |"
    $Dec = ("+", ("=" * ($Str.Length - 2)), "+" -Join "")
    write-Host $Dec`n$Str`n$Dec
}