using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Tools.IO;


namespace Tools.Converter
{
    public class ConvertingRoutines
    {

        /// <summary>
        /// Creates a LinkedList containing the rows for a CSV file.
        /// </summary>
        /// <param name="values">The values for the CSV file.</param>
        /// <param name="separator">How do we separate values in the same row?</param>
        /// <param name="separator2">Marks the end of the row.</param>
        /// <returns>List of strings that can be saved as a CSV file.</returns>
        public static LinkedList<string> List2Strings<S, T>(ICollection<S> values, string separator, string separator2) where S : ICollection<T>
        {
            LinkedList<string> res = new LinkedList<string>();
            StringBuilder str = new StringBuilder();
            foreach (S row in values)
            {
                str.Clear();
                foreach (T value in row)
                    str.Append(value).Append(separator);
                res.AddLast(str.ToString(0, str.Length - separator.Length) + separator2);
            }
            return res;
        }

        /// <summary>
        /// Creates a LinkedList containing the rows for a CSV file.
        /// </summary>
        /// <param name="values">The values for the CSV file.</param>
        /// <param name="separator">How do we separate values in the same row?</param>
        /// <returns>List of strings that can be saved as a CSV file.</returns>
        public static LinkedList<string> List2Strings<S, T>(ICollection<S> values, string separator) where S : ICollection<T>
        {
            return List2Strings<S, T>(values, separator, "");
        }
    }
}
