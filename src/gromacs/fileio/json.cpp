/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#include "json.h"
#include <algorithm>
#include <cctype>

namespace json
{

std::string pull_front_key_from_string(std::string &s)
{
    std::string::iterator key_start = find(s.begin(), s.end(), '\"');
    if (key_start == s.end())
    {
        throw std::out_of_range("JSON string parsing : Cannot find opening quotes for key in \'"+ s + "\'.");
    }
    std::string::iterator key_end  =  find(key_start+1, s.end(), '\"');
    if (key_end == s.end())
    {
        throw std::out_of_range("JSON string parsing : Cannot find closing quotes for key in \'"+ s + "\'.");
    }
    std::string key(key_start+1, key_end);
    s.erase(s.begin(), key_end+1);
    return key;
};

std::string
pull_matching_bracket_string_from_string(std::string &s, const char open_bracket, const char close_bracket)
{
    int                   bracket_balance           = 1;
    std::string::iterator opening_bracket_position  = find(s.begin(), s.end(), open_bracket);
    std::string::iterator current_position          = opening_bracket_position;
    std::string::iterator next_close_position;
    std::string::iterator next_open_position;
    while (bracket_balance != 0)
    {
        next_close_position = find(current_position, s.end(), close_bracket);
        next_open_position  = find(current_position, s.end(), open_bracket);
        if (next_close_position == s.end())
        {
            throw std::out_of_range("Error parsing brackets: Missing closing bracket.");
        }
        if (next_open_position < next_close_position)
        {
            ++bracket_balance;
            current_position = next_open_position;
        }
        else
        {
            --bracket_balance;
            current_position = next_close_position;
        }

    }

    std::string result(opening_bracket_position, current_position+1);
    s.erase(opening_bracket_position, current_position+1);
    return result;
}

std::string
Object::operator[](const std::string &key)
{
    return value_[key]->write();
}

std::string
Object::at(const std::string &key)
{
    return value_.at(key)->write();
}


std::string
pull_array_string_from_string(std::string &s)
{
    return pull_matching_bracket_string_from_string(s, '[', ']');
};

std::string
pull_object_string_from_string(std::string &s)
{
    return pull_matching_bracket_string_from_string(s, '{', '}');
};

std::string
pull_plain_value_string_from_string(std::string &s)
{
    std::string::iterator next_comma_or_end = find(s.begin(), s.end(), ',');
    std::string           result(s.begin(), next_comma_or_end);
    if (next_comma_or_end != s.end())
    {
        s.erase(s.begin(), next_comma_or_end);
    }
    else
    {
        s.clear();
    }
    return result;
}

bool
starts_array(const std::string &s)
{
    return s.front() == '[';
}

bool
starts_object(const std::string &s)
{
    return s.front() == '{';
}

void
strip_enclosing_brackets(std::string &s, const char opening_bracket, const char closing_bracket)
{
    s.erase(s.begin(), s.begin()+s.find(opening_bracket)+1);
    s.erase(s.begin()+s.rfind(closing_bracket), s.end());
}

std::string pull_front_value_from_string(std::string &s)
{
    std::string result;
    if (starts_array(s))
    {
        return pull_array_string_from_string(s);
    }
    if (starts_object(s))
    {
        return pull_object_string_from_string(s);
    }
    return pull_plain_value_string_from_string(s);
};

bool
starts_bool(std::string s)
{
    return (s.compare(0, 4, "true") == 0) || (s.compare(0, 5, "false") == 0);
}

bool starts_number(std::string s)
{
    return std::isdigit(s[0]) || (s.compare(0, 1, ".") == 0);
}

bool starts_null(std::string s)
{
    return s.compare(0, 4, "null") == 0;
}

std::unique_ptr<Entry> create_value_object_from_string(std::string &s)
{
    if (starts_array(s))
    {
        return std::unique_ptr<Entry>(new Array(s));
    }
    if (starts_object(s))
    {
        return std::unique_ptr<Entry>(new Object(s));
    }
    if (starts_bool(s))
    {
        return std::unique_ptr<Entry>(new Bool(s));
    }
    if (starts_number(s))
    {
        return std::unique_ptr<Entry>(new Number(s));
    }
    if (starts_null(s))
    {
        return std::unique_ptr<Entry>(new Null(s));
    }
    return std::unique_ptr<Entry>(new String(s));
}

Object::Object(std::string &s)
{
    // remove all whitespace
    s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );
    strip_enclosing_brackets(s, '{', '}');
    std::string key;
    while (s.find("\"") != std::string::npos)
    {
        std::string key          = pull_front_key_from_string(s);
        s.erase(s.begin(), s.begin()+s.find(":")+1);
        std::string value_string = pull_front_value_from_string(s);
        s.erase(s.begin(), s.begin()+s.find(",")+1);
        value_[key] = create_value_object_from_string(value_string);
    }

};

std::string Object::write()
{
    std::string result("{");
    for (auto &item : value_)
    {
        result += std::string("\"") + item.first +   std::string("\"") + std::string(":") + item.second->write() + std::string(",");
    }
    return result + std::string("}");
};

Array::Array(std::string &s)
{
    std::string value_string;
    strip_enclosing_brackets(s, '[', ']');
    while (s.find(",") != std::string::npos)
    {
        value_string = pull_front_value_from_string(s);
        value_.push_back(create_value_object_from_string(value_string));
        s.erase(s.begin(), s.begin()+s.find(",")+1);
    }
};

std::string Array::write()
{
    std::string result("[");
    for (auto &item : value_)
    {
        result += item->write()+std::string(",");
    }
    return result + std::string("]");
};

Number::Number(std::string &s)
{
    char ** end = nullptr;
    value_ = std::strtof(s.c_str(), end);
}
std::string Number::write()
{
    return std::string(std::to_string(value_));
}

String::String(std::string &s)
{
    strip_enclosing_brackets(s, '\"', '\"');
    value_ = s;
};

std::string String::write()
{
    return std::string(value_);
}
;

Bool::Bool(std::string &s)
{
    if (s.find("true") != std::string::npos)
    {
        value_ = true;
    }
    else
    {
        value_ = false;
    }
};

std::string Bool::write()
{
    return value_ ? std::string("true") : std::string("false");
};

Null::Null(std::string & /*s*/)
{
}
std::string Null::write()
{
    return std::string("null");
}

}
