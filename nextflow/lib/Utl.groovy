import groovy.json.JsonOutput

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

    public static delSuffix(ch) {
        return ch.map({ [it[0] + ["suffix": null]] + it[1..-1] })
    }

    public static getId(ch, prepend = false) {
        if (prepend) {
            return ch.map({ [it[0].id] + it })
        } else {
            return ch.map({ [it[0].id] + [it[1]] })
        }
    }

    public static delId(ch) {
        return ch.map({ it[1..-1] })
    }

    public static addSuffix(ch, suffix) {
        return ch.map({
            def val = suffix instanceof Closure ? suffix(it[0]) : suffix
            [it[0] + ["suffix": val]] + it[1..-1]
        })
    }

    /**  Change the meta map (first element) of a given channel
    /** The values of `override` can be literals or closures.    **/
    /**  In the latter case, they must be a closure that takes the old meta as input    **/
    /**     and returns a literal value
    **/
    public static modifyMeta(ch, Map override) {
        def changeVals = {
            def processed = [:]
            override.each({ k, v ->
                if (v instanceof Closure) {
                    processed.put(k, v(it[0]))
                } else {
                    processed.put(k, v)
                }
            })
           [it[0] + processed] + it[1..-1]
        }
        return ch.map(changeVals)
    }

    public static mapToJson(Map m) {
        return JsonOutput.prettyPrint(JsonOutput.toJson(m))
    }

    // Remove cli args or flags specified in args_a that have been specified in args_b
    public static overrideArgs(args_a, args_b) {
        def to_check = args_b.collect({ it.split()[0] })
        def final_args = args_a.findAll({ !to_check.contains(it.split()[0]) })
        return final_args.join(" ")
    }

    // Convert a list of maps `lst` into a string csv representation
    public static listOfMaps2Csv(lst) {
        def keys = [lst[0].keySet().join(",")]
        def joined = lst.collect { it.values().join(",") }
        return (keys + joined).join("\n")
    }

}
