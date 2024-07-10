
using module ./CMake.psm1
using module ./Git.psm1
using module ./Utilities.psm1

class GenericInstaller {
    hidden [string] $Name
    hidden [string] $Repository
    hidden [string] $Tag
    hidden [string] $BuildType
    hidden [int] $ParallelJobs = 1
    hidden [string[]] $Options = @()

    hidden [string] $RepositoryDir
    hidden [string] $BuildDir
    hidden [string] $InstallDir

    hidden [bool] $DoInstall = $false  

    hidden [Git] $Git
    hidden [CMake] $CMake
    
    GenericInstaller([hashtable] $Properties) {
        [GenericInstaller]::ValidateProperties($Properties)

        foreach ($Property in $Properties.Keys) {
            $this.$Property = $Properties.$Property
        }

        $this.ProcessDirectories()

        $this.Git = [Git]::new(@{
                Name          = $this.Name
                Repository    = $this.Repository
                Tag           = $this.Tag
                RepositoryDir = $this.RepositoryDir
            })

        $this.CMake = [CMake]::new(@{
                Name          = $this.Name
                RepositoryDir = $this.RepositoryDir
                BuildDir      = $this.BuildDir
                InstallDir    = $this.InstallDir
                BuildType     = $this.BuildType
                ParallelJobs  = $this.ParallelJobs
                Options       = $this.Options
            })
    }

    hidden static [void] ValidateProperties([hashtable] $Properties) {

        [string[]] $requiredProperties = @(
            'Name', 
            'Repository', 
            'Tag',
            'BuildType'
        )

        [string[]] $optionalProperties = @(
            'ParallelJobs',
            'Options'
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

    hidden [void] ProcessDirectories() {
        $this.InstallDir = (Join-Path $Global:InstallDir $this.Name)
        if (-not(Test-Path -Path $this.InstallDir)) {
            New-Item $this.InstallDir -Type Directory
            $this.RepositoryDir = (Join-Path $Global:RepositoryDir $this.Name)
            New-Item $this.RepositoryDir -Type Directory
            $this.BuildDir = (Join-Path $Global:BuildDir $this.Name)
            New-Item $this.BuildDir -Type Directory  
            $this.DoInstall = $true
        }
    }

    [string] Directory() { return $this.InstallDir }

    [string] StaticLibraryDirectory() { return (Join-Path  $this.InstallDir 'lib') }

    [string] IncludeDirectory() { return (Join-Path  $this.InstallDir 'include') }

    [void] Run() {
        if ($this.DoInstall) {
            Utilities\PrintDecoratedMessage($this.Name, "Start")
            $this.Git.Run()
            $this.CMake.Run()
            Utilities\PrintDecoratedMessage($this.Name, "End")
        }
        else {
            Write-Host $this.Name is already installed.
        }
    }
}
