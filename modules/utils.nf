/* Shared helper functions */

// Same location every run, so it never depends on the run name.
def runSettingsPath(launchDir) {
    return launchDir.resolve('.nextflow').resolve('run.settings')
}

// True when a SQLAlchemy URL embeds a password (user:pass@host).
def urlHasUserInfo(String url) {
    return url ? url.matches('^[^:]+://[^/@:]+:[^/@]*@.*$') : false
}

// Write the run's settings file, owner read/write only.
def writeRunSettings(path, String dbUrl, String apiKey) {
    java.nio.file.Files.createDirectories(path.parent)
    def ownerOnly = java.nio.file.attribute.PosixFilePermissions.fromString('rw-------')
    java.nio.file.Files.deleteIfExists(path)
    try {
        java.nio.file.Files.createFile(path, java.nio.file.attribute.PosixFilePermissions.asFileAttribute(ownerOnly))
    } catch (UnsupportedOperationException _e) {
        java.nio.file.Files.createFile(path)                       // non-POSIX filesystem
    }
    def text = "DB_URL=${dbUrl ?: ''}\nAPI_KEY=${apiKey ?: ''}\n"
    java.nio.file.Files.write(path, text.getBytes('UTF-8'), java.nio.file.StandardOpenOption.TRUNCATE_EXISTING)
}

def deleteRunSettings(path) {
    java.nio.file.Files.deleteIfExists(path)
}

// Tasks only get this path, single-quoted for bash.
def runSettingsArg(launchDir) {
    return "--settings-file '" + runSettingsPath(launchDir).toString().replace("'", "'\"'\"'") + "'"
}
