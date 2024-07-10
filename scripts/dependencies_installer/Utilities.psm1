Function PrintDecoratedMessage([string] $Name, [string] $Message) {
    $Str = "| $Name Installer: $Message |"
    $Dec = ("+", ("=" * ($str.Length - 2)), "+" -Join "")
    write-Host $Dec`n$Str`n$Dec
}