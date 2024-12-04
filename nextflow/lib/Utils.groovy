class Utils {
    public static String getName(Integer module_number, Map meta, String suffix = null, String ext = null) {
        def suf = meta.suffix ? meta.suffix : suffix
        if (!ext && suffix) {
            return "${module_number}-${meta.filename}-${suf}"
        } else if (ext && suf) {
            return "${module_number}-${meta.filename}-${suf}.${ext}"
        } else if (!suf && !ext) {
            return "${module_number}-${meta.filename}"
        } else {
            return "${module_number}-${meta.filename}.${ext}"
        }
    }

    // public static String saveFn(x, String patterns = /.*\.log/) {
    //     x ==~ patterns ||
    // }
}
