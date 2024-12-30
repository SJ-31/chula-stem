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

    public static joinById(first, to_join) {
        def prepend_id = first.map({ [it[0].id] + it })
        to_join.each {
            prepend_id = prepend_id.join(it.map({ [it[0].id] + it[1..-1] }))
        }
        return prepend_id.map({ it[1..-1] })
    }

    public static getId(x) {
        return [x[0].id] + [x[1]]
    }
}
