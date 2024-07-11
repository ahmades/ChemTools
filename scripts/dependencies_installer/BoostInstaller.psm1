using module ./Utilities.psm1

class BoostInstaller {
    hidden [string] $Tag
    hidden [string] $BuildType
    hidden [int] $ParallelJobs = 1

    hidden [string] $Name = 'boost'
    hidden [string] $BuildDir = (Join-Path $Global:BuildDir $this.Name)
    hidden [string] $InstallDir = (Join-Path $Global:InstallDir $this.Name)
    hidden [string] $ExtractDir
    hidden [string] $ArchivePath
    hidden [string] $Stem

    hidden [bool] $DoInstall = $false

    # Default constructor
    BoostInstaller() { $this.Initialise(@{}) }
    
    BoostInstaller([hashtable] $Properties) {
        [BoostInstaller]::ValidateProperties($Properties)
        $this.Initialise($Properties)
        $this.ProcessDirectories()
    }

    hidden static [void] ValidateProperties([hashtable] $Properties) {

        [string[]] $requiredProperties = @(
            'Tag',
            'BuildType'
        )

        [string[]] $optionalProperties = @(
            'ParallelJobs'
        )

        $allProperties = $requiredProperties + $optionalProperties

        foreach ($property in $requiredProperties) {
            if (-not $Properties.ContainsKey($property)) {
                throw "Not all required properties are provided"
            }
        }

        foreach ($property in $Properties.Keys) {
            if ($allProperties -notcontains $property) {
                throw "Invalid properties are provided"
            }
        }
    }

    # Initialiser method
    hidden [void] Initialise([hashtable]$Properties) {
        foreach ($Property in $Properties.Keys) {
            $this.$Property = $Properties.$Property
        }
    }

    hidden [void] ProcessDirectories() {
        if (-not(Test-Path -Path $this.InstallDir)) {
            New-Item $this.InstallDir -Type Directory
            if (Test-Path -Path $this.BuildDir) {
                Remove-Item $this.BuildDir -Recurse -Force
            }
            New-Item $this.BuildDir -Type Directory
            $this.DoInstall = $true
        }
    }

    hidden [void] Downlaod() {
        $VersionUnderScore = $this.Tag.Replace('.', '_')
        $this.Stem = $this.Name + '_' + $VersionUnderScore
        $FileName = $this.Stem + ($Global:IsLinux ? '.tar.gz' : '.zip')
        $this.ArchivePath = (Join-Path $Global:DownloadDir $FileName)

        if (-not(Test-Path -Path $this.ArchivePath)) {
            $URL = (@('https://boostorg.jfrog.io/artifactory/main/release'; $this.Tag; 'source'; $FileName) -Join '/')
            $WebClient = [System.Net.WebClient]::new()
            $WebClient.DownloadFile($URL, $this.ArchivePath)
            $ExitCode = $LASTEXITCODE
            if ( $ExitCode -ne 0 ) {
                throw ('Failed to download boost.')
            }
        }
    }

    hidden [void] Extract() {
        $this.ExtractDir = (Join-Path $Global:ExtractDir $this.Stem)
        if (-not(Test-Path -Path $this.ExtractDir)) {
            if ($Global:IsLinux) {
                New-Item $this.ExtractDir -Type Directory
                tar -xvf $this.ArchivePath -C $Global:ExtractDir
            }
            else {
                Add-Type -Assembly "System.IO.Compression.Filesystem"
                [System.IO.Compression.ZipFile]::ExtractToDirectory($this.ArchivePath, $Global:ExtractDir)
                #Expand-Archive -Force  $this.ArchivePath -DestinationPath $Global:ExtractDir
            }
            $ExitCode = $LASTEXITCODE
            if ( $ExitCode -ne 0 ) {
                throw ('Failed to extract boost.')
            }
        }
    }

    hidden [void] Bootstrap() {
        Push-Location $this.ExtractDir
        Start-Process -FilePath ('./bootstrap' + ($Global:IsLinux ? '.sh' : '.bat')) -Wait -NoNewWindow
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -ne 0 ) {
            throw ('Failed to boostrap boost.')
        }
        Pop-Location
    }

    hidden [void] Install() {
        Push-Location $this.ExtractDir
        $BOOSTBuildArgumentList = `
            ' -j' + $this.ParallelJobs `
            + ' --build-dir=' + $this.BuildDir `
            + ' --prefix=' + $this.InstallDir `
            + ' variant=release runtime-link=static link=static address-model=64 install'
        Start-Process -FilePath './b2' -ArgumentList $BOOSTBuildArgumentList -Wait -NoNewWindow
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -ne 0 ) {
            throw ('Failed to install boost.')
        }
        Pop-Location
    }

    [string] Directory() { return $this.InstallDir }

    [string] StaticLibraryDirectory() { return (Join-Path  $this.InstallDir 'lib') }

    [string] IncludeDirectory() { 
        $VersionUnderScore = $this.Tag.Replace('.', '_')
        $IndexOfLastUnderscore = $VersionUnderScore.LastIndexOf("_")
        $VersionUnderScoreMajorMinor = $VersionUnderScore.Substring(0, $IndexOfLastUnderscore)
        return (Join-Path  $this.InstallDir 'include' ($this.Name + '-' + $VersionUnderScoreMajorMinor)) 
    }

    [void] Run() {
        if ($this.DoInstall) {
            Utilities\PrintDecoratedMessage $this.Name "Start"
            $this.Downlaod()
            $this.Extract()
            $this.Bootstrap()
            $this.Install()
            Utilities\PrintDecoratedMessage $this.Name "End"
        }
        else {
            Write-Host $this.Name is already installed.
        }
    }
}
