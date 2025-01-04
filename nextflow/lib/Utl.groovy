class Utl {
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

  /** Sequentially join all channels in list `to_join` with channel `first`     */
    /**
    /** All channels involved are expected to emit items that start with a meta map     */
    /**     Only the meta map of `first` is retained by default
     * @param on the key(s) of the meta map to join channel elements with
     */
    public static joinFirst(first, to_join, on = ["id"], keep_meta = false) {
        def n_keys = on.size()
        def by = (0..(n_keys - 1)).collect()
        // If on == ["id"], by = [0]; if on == ["id", "bar"], by = [0, 1]

        def getVals = { meta -> on.collect { meta[it] } }
        // closure to extract all keys specified by `on`
        //  from the meta map, collecting them into a list

        def key_prepended = first.map({ getVals(it[0]) + it })

        def n = keep_meta ? 0 : 1
        to_join.each {
            key_prepended = key_prepended.join(it.map({ getVals(it[0]) + it[n..-1] }),
                                           by: by)
        }
        return key_prepended.map({ it[n_keys..-1] })
    }

    public static getId(item) {
        return [it[0].id + [it[1]]]
    }

    public static delId(item) {
        return it[1..-1]
    }

    public static prependId(item) {
        return [item[0].id] + it
    }

    public static delSuffix(ch) {
        return ch.map({ [it[0] + ["suffix": null]] + it[1..-1] })
    }

    public static addSuffix(ch, suffix) {
        return ch.map({ [it[0] + ["suffix": suffix]] + it[1..-1] })
    }

    public static getId(ch) {
        // Gets id while removing meta map
        ch.map({ [it[0].id] + it[1..-1] })
    }

    public static modifyMeta(ch, Map meta) {
        return ch.map({[it[0] + meta] + it[1..-1] })
    }

}
